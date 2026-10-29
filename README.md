# FLASHApp

FLASHApp for visualizing FLASHDeconv's results \
This app is based on [OpenMS streamlit template project](https://github.com/OpenMS/streamlit-template).

![Overview](https://github.com/user-attachments/assets/35fe2c24-7dbc-40cd-b8b5-b7504768ade1)

run locally:

`streamlit run app.py local`

### Visualizations

The viewers are built from [openms-insight](https://github.com/t0mdavid-m/openms-insight)
components (tables, heatmaps, spectra, 3D plots and the proteoform sequence view).
Until the openms-insight release with these components is on PyPI, install it
from source (the Docker and Windows builds do this for you):

```
git clone -b claude/project-thread-hxjbbs https://github.com/t0mdavid-m/openms-insight.git
cd openms-insight/js-component && npm ci && npm run build
mkdir -p ../openms_insight/js-component && cp -r dist ../openms_insight/js-component/
cd .. && pip install --no-deps .
```

## Build

Set streamlit to production in `./.streamlit/config.toml` (`developmentMode = false`)
before building any version of FLASHApp.

### Docker

First you need to build an image locally.

build image with: `docker build -f Dockerfile --no-cache -t flashapp:latest --build-arg GITHUB_TOKEN=<your-github-token> .`

You should see a successful output, but you can check if an image is built with:

`docker image ls`

After it has been built you can run the image with:

`docker run -p 8501:8501 flashapp:latest`

Navigate to `http://localhost:8501` in your browser.

## Legal pages (Impressum, Privacy Policy, Terms of Use)

Every page shows **Impressum**, **Privacy Policy** and **Terms of Use** links in the
sidebar footer, and the GDPR consent banner links to the privacy policy. By default
these point to the official OpenMS pages (`https://openms.de/impressum`, `/privacy`,
`/terms`). To override them — for example when self-hosting or deploying FLASHApp
under a different operator — set `legal_links` in `settings.json`:

    "legal_links": {
        "impressum": "https://your-domain.example/impressum",
        "privacy": "https://your-domain.example/privacy",
        "terms": "https://your-domain.example/terms"
    }

Any link you omit falls back to its OpenMS default. The `privacy` URL is reused for the
consent banner's privacy-policy link, so consent and policy stay in sync.
