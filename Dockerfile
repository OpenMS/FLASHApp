# This Dockerfile builds OpenMS, the TOPP tools, pyOpenMS and thidparty tools.
# It also adds a basic streamlit server that serves a pyOpenMS-based app.
# hints:
# build image and give it a name (here: streamlitapp) with: docker build --no-cache -t streamlitapp:latest --build-arg GITHUB_TOKEN=<your-github-token> . 2>&1 | tee build.log
# check if image was build: docker image ls
# run container: docker run -p 8501:8501 streamlitappsimple:latest
# debug container after build (comment out ENTRYPOINT) and run container with interactive /bin/bash shell
# prune unused images/etc. to free disc space (e.g. might be needed on gitpod). Use with care.: docker system prune --all --force


FROM ubuntu:24.04 AS setup-build-system
WORKDIR /

ARG OPENMS_REPO=https://github.com/t0mdavid-m/OpenMS.git
ARG OPENMS_BRANCH=FVdeploy
ARG PORT=8501
# Streamlit app Gihub user name (to download artifact from).
ARG GITHUB_USER=OpenMS
# Streamlit app Gihub repository name (to download artifact from).
ARG GITHUB_REPO=FLASHApp
# Name of the zip file containing the windows executable
ARG ASSET_NAME=OpenMS-App.zip


USER root

# Install required Ubuntu packages.
RUN apt-get -y update
RUN apt-get install -y --no-install-recommends --no-install-suggests g++ autoconf automake patch libtool make git gpg wget ca-certificates curl jq libgtk2.0-dev openjdk-8-jdk cron cmake
RUN update-ca-certificates
RUN apt-get install -y --no-install-recommends --no-install-suggests libsvm-dev libeigen3-dev coinor-libcbc-dev libglpk-dev libzip-dev zlib1g-dev libxerces-c-dev libbz2-dev libomp-dev libhdf5-dev
RUN apt-get install -y --no-install-recommends --no-install-suggests libboost-date-time-dev \
                                                                     libboost-iostreams-dev \
                                                                     libboost-regex-dev \
                                                                     libboost-math-dev \
                                                                     libboost-random-dev
RUN apt-get install -y --no-install-recommends --no-install-suggests qt6-base-dev libqt6svg6-dev libqt6opengl6-dev libqt6openglwidgets6 libgl-dev

# Install Github CLI
RUN (type -p wget >/dev/null || (apt-get update && apt-get install wget -y)) \
	&& mkdir -p -m 755 /etc/apt/keyrings \
	&& wget -qO- https://cli.github.com/packages/githubcli-archive-keyring.gpg | tee /etc/apt/keyrings/githubcli-archive-keyring.gpg > /dev/null \
	&& chmod go+r /etc/apt/keyrings/githubcli-archive-keyring.gpg \
	&& echo "deb [arch=$(dpkg --print-architecture) signed-by=/etc/apt/keyrings/githubcli-archive-keyring.gpg] https://cli.github.com/packages stable main" | tee /etc/apt/sources.list.d/github-cli.list > /dev/null \
	&& apt-get update \
	&& apt-get install gh -y

# Download and install miniforge.
ENV PATH="/root/miniforge3/bin:${PATH}"
RUN wget -q \
    https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh \
    && bash Miniforge3-Linux-x86_64.sh -b \
    && rm -f Miniforge3-Linux-x86_64.sh
RUN mamba --version

# Make /root traversable so the entrypoint can `source
# /root/miniforge3/bin/activate ...` when the container runs as a non-root
# user (apptainer/singularity maps the host UID into the container; the
# default ubuntu /root is 0700 which would block path traversal). +x only,
# not +r, so the directory listing remains private.
RUN chmod o+x /root

# Setup mamba environment.
RUN mamba create -n streamlit-env python=3.11
RUN echo "mamba activate streamlit-env" >> ~/.bashrc
SHELL ["/bin/bash", "--rcfile", "~/.bashrc"]
SHELL ["mamba", "run", "-n", "streamlit-env", "/bin/bash", "-c"]

RUN pip install --upgrade pip && python -m pip install -U setuptools nose 'Cython>=3.1' 'autowrap==0.24' pandas 'numpy>=2.0' pytest

# Clone OpenMS branch and the associcated contrib+thirdparties+pyOpenMS-doc submodules.
RUN git clone --recursive --depth=1 -b ${OPENMS_BRANCH} --single-branch ${OPENMS_REPO} && cd /OpenMS

# Pull Linux compatible third-party dependencies and store them in directory thirdparty.
WORKDIR /OpenMS
RUN mkdir /thirdparty && \
    git submodule update --init THIRDPARTY && \
    cp -r THIRDPARTY/All/* /thirdparty && \
    cp -r THIRDPARTY/Linux/x86_64/* /thirdparty && \
    chmod -R +x /thirdparty
ENV PATH="/thirdparty/LuciPHOr2:/thirdparty/MSGFPlus:/thirdparty/ThermoRawFileParser:/thirdparty/Comet:/thirdparty/Percolator:/thirdparty/Sage:${PATH}"

# Build the OpenMS-Insight package (Python + Vue bundle) from the migration branch.
# Insight's Vue dist is gitignored and it has no pip build hook, so build the bundle
# here and sync it into the package tree; the compile-openms stage pip-installs it.
FROM node:21 AS insight-build
ARG INSIGHT_REPO=https://github.com/t0mdavid-m/openms-insight.git
ARG INSIGHT_BRANCH=claude/kind-heisenberg-u6dVm
ADD https://api.github.com/repos/t0mdavid-m/openms-insight/git/refs/heads/${INSIGHT_BRANCH} insight-ref.json
RUN git clone -b ${INSIGHT_BRANCH} --single-branch ${INSIGHT_REPO} /openms-insight
WORKDIR /openms-insight/js-component
RUN npm install && npm run build
RUN mkdir -p /openms-insight/openms_insight/js-component \
    && rm -rf /openms-insight/openms_insight/js-component/dist \
    && cp -r /openms-insight/js-component/dist /openms-insight/openms_insight/js-component/dist \
    && rm -rf /openms-insight/js-component/node_modules

# Build OpenMS and pyOpenMS.
FROM setup-build-system AS compile-openms
WORKDIR /

# Set up build directory.
RUN mkdir /openms-build
WORKDIR /openms-build

# Configure (two-pass: first without miniforge so cmake finds system C++ libs,
# then reconfigure with miniforge for pyopenms).
SHELL ["/bin/bash", "-c"]
RUN cmake -DCMAKE_BUILD_TYPE='Release' -DCMAKE_PREFIX_PATH='/OpenMS/contrib-build/;/usr/;/usr/local' -DCMAKE_IGNORE_PREFIX_PATH=/root/miniforge3 -DHAS_XSERVER=OFF -DBOOST_USE_STATIC=OFF -DOPENMP=ON ../OpenMS
SHELL ["mamba", "run", "-n", "streamlit-env", "/bin/bash", "-c"]
RUN cmake -DPYOPENMS=ON -DPY_MEMLEAK_DISABLE=On -DCMAKE_IGNORE_PREFIX_PATH=/root/miniforge3 .

# Build TOPP tools and clean up.
RUN make -j4 TOPP
RUN rm -rf src doc

# Build pyOpenMS wheels and install via pip.
RUN make -j4 pyopenms
WORKDIR /openms-build/pyOpenMS
RUN pip install dist/*.whl

# Install other dependencies (excluding pyopenms)
COPY requirements.txt ./requirements.txt 
RUN grep -Ev '^pyopenms([=<>!~].*)?$' requirements.txt > requirements_cleaned.txt && mv requirements_cleaned.txt requirements.txt
# OpenMS-Insight: install from the migration branch built in the insight-build stage
# (with its Vue bundle), before requirements so the pin resolves from source not PyPI.
COPY --from=insight-build /openms-insight /tmp/openms-insight
RUN pip install /tmp/openms-insight && rm -rf /tmp/openms-insight
RUN pip install -r requirements.txt

WORKDIR /
RUN mkdir openms

# Copy TOPP tools bin directory, add to PATH.
RUN cp -r openms-build/bin /openms/bin
ENV PATH="/openms/bin/:${PATH}"

# Copy TOPP tools bin directory, add to PATH.
RUN cp -r openms-build/lib /openms/lib
ENV LD_LIBRARY_PATH="/openms/lib/:/root/miniforge3/envs/streamlit-env/lib:${LD_LIBRARY_PATH}"

# Copy share folder, add to PATH, remove source directory.
RUN cp -r OpenMS/share/OpenMS /openms/share
RUN rm -rf OpenMS
ENV OPENMS_DATA_PATH="/openms/share/"

# Remove build directory.
RUN rm -rf openms-build

# Prepare and run streamlit app.
# (The legacy local Vue component build stage was removed in the OpenMS-Insight
# migration -- Insight ships its own Vue bundle via the openms-insight package.)
FROM compile-openms AS run-app

# Install Redis server for job queue and nginx for load balancing.
# Redis data lives under $RUNTIME_DIR at runtime (see entrypoint.sh) so no
# /var/lib/redis setup is needed - that path is not writable under Apptainer.
RUN apt-get update && apt-get install -y --no-install-recommends redis-server nginx \
    && rm -rf /var/lib/apt/lists/*

# Create Redis data directory. Default 0755 root-owned is enough: the docker
# entrypoint runs as root (can write regardless of mode), and the apptainer
# entrypoint relocates Redis state to /tmp/openms-runtime-* so this dir is
# never written under apptainer.
RUN mkdir -p /var/lib/redis

# Pre-create bind-mount targets so apptainer/singularity has a real attach
# point. Docker auto-creates missing `-v` targets, but singularity uses a
# read-only underlay and silently ignores `:rw` when the target isn't a
# real directory in the SIF — writes then fail with EROFS even though the
# host bind path is writable. Pre-creating these directories costs one
# inode each and changes nothing in docker mode (the user's volume mount
# shadows them).
RUN mkdir -p /workspaces-streamlit-template /mounted-data

# Create workdir and copy over all streamlit related files/folders.

# note: specifying folder with slash as suffix and repeating the folder name seems important to preserve directory structure
WORKDIR /app

COPY .streamlit/ /app/.streamlit
COPY assets/ /app/assets
COPY static/ /app/static
COPY clean-up-workspaces.py /app/clean-up-workspaces.py
COPY content/ /app/content
COPY example-data/ /app/example-data
COPY gdpr_consent/ /app/gdpr_consent
COPY hooks/ /app/hooks
COPY src/ /app/src
COPY app.py /app/app.py
COPY settings.json /app/settings.json
COPY default-parameters.json /app/default-parameters.json
COPY presets.json /app/presets.json

# add cron job to the crontab
RUN echo "0 3 * * * /root/miniforge3/envs/streamlit-env/bin/python /app/clean-up-workspaces.py >> /app/clean-up-workspaces.log 2>&1" | crontab -

# Set default worker count (can be overridden via environment variable)
ENV RQ_WORKER_COUNT=1
ENV REDIS_URL=redis://localhost:6379/0

# Number of Streamlit server instances for load balancing (default: 1 = no load balancer)
# Set to >1 to enable nginx load balancer with multiple Streamlit instances
ENV STREAMLIT_SERVER_COUNT=1

# Install the apptainer-compatible entrypoint that starts cron (when the root
# FS is writable), Redis, RQ workers, optional nginx load balancer, and the
# Streamlit server. The script falls back to /tmp paths under apptainer.
COPY docker/entrypoint.sh /app/entrypoint.sh
RUN chmod +x /app/entrypoint.sh

# Patch Analytics
RUN mamba run -n streamlit-env python hooks/hook-analytics.py

# Set Online Deployment
RUN jq '.online_deployment = true' settings.json > tmp.json && mv tmp.json settings.json

# Download latest OpenMS App executable as a ZIP file
# Re-declare ARGs needed in this stage (ARGs don't persist across FROM)
ARG GITHUB_TOKEN
ARG GITHUB_USER=OpenMS
ARG GITHUB_REPO=FLASHApp
# RELEASE_TAG pins the download to the release being published (set by the
# build-and-test workflow on release events). When empty we fall back to the
# latest release, preserving the previous behavior for develop/manual builds.
ARG RELEASE_TAG
RUN if [ -n "$GITHUB_TOKEN" ]; then \
        echo "Downloading release asset..."; \
        if [ -n "$RELEASE_TAG" ]; then \
            GH_TOKEN="$GITHUB_TOKEN" gh release download "$RELEASE_TAG" -R ${GITHUB_USER}/${GITHUB_REPO} -p "OpenMS-App.zip" -D /app; \
        else \
            GH_TOKEN="$GITHUB_TOKEN" gh release download -R ${GITHUB_USER}/${GITHUB_REPO} -p "OpenMS-App.zip" -D /app; \
        fi; \
    else \
        echo "No token, skipping download."; \
    fi


# Run app as container entrypoint.
EXPOSE $PORT
ENTRYPOINT ["/app/entrypoint.sh"]
