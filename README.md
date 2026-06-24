# FLASHApp

FLASHApp for visualizing FLASHDeconv's results \
This app is based on [OpenMS streamlit template project](https://github.com/OpenMS/streamlit-template).

![Overview](https://github.com/user-attachments/assets/35fe2c24-7dbc-40cd-b8b5-b7504768ade1)

run locally:

`streamlit run app.py local`

### Working with submodules

This project uses a git submodule to integrate openms-streamlit-vue-component.

When checking out this repository you will need to run a few extra commands:

`git submodule init`

`git submodule update`

If you would like to update the submodule to the latest commit, use the following:

`git submodule update --remote`

See [https://git-scm.com/book/en/v2/Git-Tools-Submodules](https://git-scm.com/book/en/v2/Git-Tools-Submodules)
for more documentation on submodules

## Build

To build FLASHApp, you first need to build the `openms-streamlit-vue-component`
and copy the build from `./openms-streamlit-vue-component/dist` to 
`./js-component/dist`.

You then should set streamlit to production in two locations:

* in `./.streamlit/config.toml` set `developmentMode` to `false`
* in `./src/components.py` set `_RELEASE` to `True`

These steps should be done before building any version of FLASHApp.

### Docker

First you need to build an image locally.

Prerequisite: `src/components.py` has `RELEASE=True` and `dist/` contains a build of the Vue
component. These should be the settings on the `main` branch.

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
