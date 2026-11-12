import streamlit as st

from src.common.common import page_setup, themed_asset, v_space
from src.common.contributors import render_contributors
from src.common.installers import render_installer_buttons

page_setup(page="main")

CLI_TOOLS_URL = "https://abibuilder.cs.uni-tuebingen.de/archive/openms/OpenMSInstaller/experimental/FVdeploy/"


def render_download_section():
    """The downloads, below the title and logo: the desktop app for each platform and
    the FLASH* command line tools."""
    # The offline downloads are pointless inside the offline app itself.
    if st.session_state.get("location") == "local":
        return
    with st.container(key="fa_downloads"):
        st.markdown(
            """
            <div class="fa-downloads-title">⬇️ Download FLASHApp</div>
            <div class="fa-downloads-text">
                Run FLASHApp on your own computer: the same analysis and results,
                on files that stay on your disk, with no upload.
            </div>
            """,
            unsafe_allow_html=True,
        )
        if not render_installer_buttons(app_name="FLASHApp"):
            st.caption("The desktop app of the latest release cannot be listed right now. "
                       "See the [releases](https://github.com/OpenMS/FLASHApp/releases).")
        st.markdown(
            f'<div class="fa-downloads-cli">Prefer scripting? Get the '
            f'<a href="{CLI_TOOLS_URL}" target="_blank">FLASH* command line tools</a> '
            f'(FLASHDeconv, FLASHTnT and more, as part of OpenMS).</div>',
            unsafe_allow_html=True,
        )
    st.markdown(
        """
        <style>
        .st-key-fa_downloads {
            border: 2px solid var(--fa-accent) !important;
            border-radius: 12px !important;
            padding: 1.25rem 1.5rem !important;
            margin: 0 0 1.5rem 0 !important;
            background: linear-gradient(135deg, rgba(128, 128, 128, 0.04) 0%, rgba(128, 128, 128, 0.10) 100%) !important;
            text-align: center !important;
        }
        .fa-downloads-title {
            color: var(--fa-accent);
            font-size: 1.5rem;
            font-weight: 700;
            margin-bottom: 0.25rem;
        }
        .fa-downloads-text { opacity: 0.85; margin-bottom: 0.75rem; }
        .fa-downloads-cli { opacity: 0.75; font-size: 0.9rem; margin-top: 0.25rem; }
        .st-key-fa_downloads [data-testid^="stBaseLinkButton"],
        .st-key-fa_downloads [data-testid^="stBaseButton"] {
            background: var(--fa-accent) !important;
            border-color: var(--fa-accent) !important;
            color: #ffffff !important;
            font-weight: 600 !important;
            min-height: 3rem !important;
        }
        .st-key-fa_downloads [data-testid^="stBaseLinkButton"] p,
        .st-key-fa_downloads [data-testid^="stBaseButton"] p {
            color: #ffffff !important;
            font-size: 1.05rem !important;
        }
        </style>
        """,
        unsafe_allow_html=True,
    )


def inject_theme_colors():
    """Define the accent and body text colors for the active light or dark theme."""
    dark = st.context.theme.type == "dark"
    st.markdown(
        f"""
        <style>
        :root {{
            --fa-accent: {"#7d8ce8" if dark else "#29379b"};
            --fa-body-text: {"#c9ccd3" if dark else "#495057"};
        }}
        </style>
        """,
        unsafe_allow_html=True,
    )


def inject_workflow_button_css():
    """Inject CSS for custom workflow button styling with responsive design."""
    st.markdown(
        """
        <style>
        /* Main workflow button styling */
        .workflow-button {
            background: linear-gradient(135deg, rgba(128, 128, 128, 0.06) 0%, rgba(128, 128, 128, 0.12) 100%);
            border: 2px solid rgba(128, 128, 128, 0.3);
            border-radius: 12px;
            padding: 2rem 1.5rem;
            text-align: center;
            cursor: pointer;
            transition: all 0.3s ease;
            height: max(280px, 20vh);
            display: flex;
            flex-direction: column;
            justify-content: center;
            align-items: center;
            text-decoration: none;
            box-shadow: 0 2px 8px rgba(0,0,0,0.1);
            margin-bottom: 1rem;
        }
        
        .workflow-button:hover {
            background: linear-gradient(135deg, #29379b 0%, #1e2a7a 100%);
            border-color: #29379b;
            color: white !important;
            transform: translateY(-4px);
            box-shadow: 0 8px 24px rgba(41, 55, 155, 0.3);
        }
        
        .workflow-button:active {
            background: linear-gradient(135deg, #1e2a7a 0%, #162159 100%);
            border-color: #1e2a7a;
            color: white !important;
            transform: translateY(-2px) scale(0.98);
            box-shadow: 0 4px 12px rgba(41, 55, 155, 0.4);
        }
        
        .workflow-button:focus {
            outline: 3px solid #29379b;
            outline-offset: 2px;
        }
        
        .workflow-button:hover .workflow-emoji {
            transform: scale(1.1);
        }
        
        .workflow-button:active .workflow-emoji {
            transform: scale(1.05);
        }
        
        .workflow-button:hover .workflow-title {
            color: white !important;
        }
        
        .workflow-button:active .workflow-title {
            color: white !important;
        }
        
        .workflow-button:hover .workflow-subtitle {
            color: #e9ecef !important;
        }
        
        .workflow-button:active .workflow-subtitle {
            color: #e9ecef !important;
        }
        
        .workflow-emoji {
            font-size: clamp(2.5rem, 4vw, 5rem);
            margin-bottom: clamp(0.75rem, 1.5vh, 1.5rem);
            transition: transform 0.3s ease;
        }
        
        .workflow-title {
            font-size: clamp(1.25rem, 2.5vw, 2rem);
            font-weight: 700;
            color: var(--fa-accent);
            margin-bottom: clamp(0.375rem, 0.75vh, 0.75rem);
            transition: color 0.3s ease;
        }
        
        .workflow-subtitle {
            font-size: clamp(0.875rem, 1.5vw, 1.25rem);
            color: #6c757d;
            font-weight: 500;
            transition: color 0.3s ease;
        }
        
        /* Enhanced download section */
        .download-section {
            background: linear-gradient(135deg, rgba(128, 128, 128, 0.06) 0%, rgba(128, 128, 128, 0.12) 100%);
            border: 2px solid rgba(128, 128, 128, 0.3);
            border-radius: 12px;
            padding: 2rem;
            margin: 2rem 0;
            text-align: center;
        }
        
        .download-section h3 {
            color: var(--fa-accent);
            margin-bottom: 1rem;
            font-size: 1.5rem;
            font-weight: 700;
        }
        
        /* Hero section styling */
        .hero-section {
            text-align: center;
            margin-bottom: 3rem;
        }
        
        .hero-title {
            font-size: 2.5rem;
            font-weight: 700;
            color: var(--fa-accent);
            margin-bottom: 1rem;
        }
        
        .hero-subtitle {
            font-size: 1.25rem;
            color: #6c757d;
            margin-bottom: 2rem;
        }
        
        /* Footer section */
        .footer-section {
            margin-top: 3rem;
            padding-top: 2rem;
            border-top: 1px solid rgba(128, 128, 128, 0.3);
        }
        
        /* Responsive design with dynamic scaling */
        
        /* Wide screens (> 1440px) - Maximum size */
        @media (min-width: 1441px) {
            .workflow-button {
                height: max(400px, 25vh);
                padding: clamp(2.5rem, 3vw, 4rem) clamp(2rem, 2.5vw, 3rem);
                max-width: 500px;
                margin: 0 auto 1.5rem auto;
            }
        }
        
        /* Desktop (1024px - 1440px) - Significantly larger */
        @media (min-width: 1024px) and (max-width: 1440px) {
            .workflow-button {
                height: max(320px, 22vh);
                padding: clamp(2rem, 2.5vw, 3rem) clamp(1.75rem, 2vw, 2.5rem);
                max-width: 450px;
                margin: 0 auto 1.25rem auto;
            }
        }
        
        /* Tablet (768px - 1023px) - Slightly larger */
        @media (min-width: 768px) and (max-width: 1023px) {
            .workflow-button {
                height: max(260px, 18vh);
                padding: clamp(1.75rem, 2vw, 2.5rem) clamp(1.5rem, 1.75vw, 2rem);
                max-width: 400px;
                margin: 0 auto 1rem auto;
            }
        }
        
        /* Mobile landscape (481px - 767px) - Moderate size */
        @media (min-width: 481px) and (max-width: 767px) {
            .workflow-button {
                height: max(240px, 16vh);
                padding: 1.5rem 1.25rem;
                margin-bottom: 1rem;
            }
            
            .hero-title {
                font-size: 2rem;
            }
            
            .hero-subtitle {
                font-size: 1.1rem;
            }
        }
        
        /* Mobile portrait (≤ 480px) - Compact size */
        @media (max-width: 480px) {
            .workflow-button {
                height: max(200px, 14vh);
                padding: 1.25rem 1rem;
                margin-bottom: 1rem;
            }
            
            .hero-title {
                font-size: 1.75rem;
            }
            
            .hero-subtitle {
                font-size: 1rem;
            }
        }
        
        /* Optimize column spacing for workflow buttons */
        .main .block-container [data-testid="column"] {
            padding-left: 0.5rem !important;
            padding-right: 0.5rem !important;
        }
        
        /* Ensure first and last columns have proper edge spacing */
        .main .block-container [data-testid="column"]:first-child {
            padding-left: 0 !important;
        }
        
        .main .block-container [data-testid="column"]:last-child {
            padding-right: 0 !important;
        }
        
        /* Reduce button margins for tighter layout */
        .workflow-button {
            margin-bottom: 0.5rem !important;
        }
        
        /* Container spacing optimization */
        .stColumn > div {
            padding-top: 0 !important;
        }
        </style>
        """,
        unsafe_allow_html=True,
    )

def render_workflow_selection():
    """Render the main workflow selection section."""
    # Hero section with title on left and FLASHApp logo on right
    st.markdown(
        """
        <div class="hero-section">
        """,
        unsafe_allow_html=True,
    )
    
    # Create columns for title and logo
    spacer1, title_col, logo_col, spacer2 = st.columns([1, 4.5, 1.5, 1])
    
    with title_col:
        st.markdown(
            """
            <h1 class="hero-title">👋 FLASHApp</h1>
            <p class="hero-subtitle">A platform for your favourite FLASH Tools!</p>
            """,
            unsafe_allow_html=True,
        )
    
    with logo_col:
        st.image(themed_asset("assets/flashapp-icon-light.svg", "assets/flashapp-icon.svg"), width=160)
    
    st.markdown("</div>", unsafe_allow_html=True)


FLOW_STEPS = [
    ("⚙️", "Describe experiment", "Pick runs, database and samples",
     "content/TopDown/TopDownConfigure.py"),
    ("🚀", "Run analysis", "Deconvolve, identify, quantify",
     "content/TopDown/TopDownRun.py"),
    ("👀", "Explore results", "Spectra, masses, proteoforms",
     "content/TopDown/TopDownViewer.py"),
    ("⬇️", "Download", "Result files of every run",
     "content/TopDown/TopDownDownload.py"),
]


def render_workflow_flow():
    """The steps of a top-down analysis, each a button to its page."""
    steps = list(FLOW_STEPS)
    if st.session_state.get("location") == "local":
        steps[-1] = ("📤", "Export", "Result files to a folder", steps[-1][3])

    st.markdown(
        '<h3 class="flow-title">⚡️ Top-down analysis</h3>',
        unsafe_allow_html=True,
    )
    widths = []
    for i in range(len(steps)):
        widths += [4, 1] if i < len(steps) - 1 else [4]
    cols = st.columns(widths, gap="small", vertical_alignment="center")

    css = []
    for i, (emoji, title, subtitle, page) in enumerate(steps):
        key = f"flow_step_{i}"
        with cols[2 * i]:
            if st.button(f"{i + 1}\\. {emoji} {title}", key=key,
                         help=subtitle, use_container_width=True):
                st.switch_page(page)
        if i < len(steps) - 1:
            with cols[2 * i + 1]:
                st.markdown('<div class="flow-arrow">➜</div>',
                            unsafe_allow_html=True)
        css.append(
            f'.st-key-{key} button::after {{ content: "{subtitle}"; }}'
        )

    st.markdown(
        "<style>"
        + """
        .flow-title {
            color: var(--fa-accent);
            text-align: center;
            font-weight: 700;
            margin: 1rem 0 1.25rem 0;
        }
        [class*="st-key-flow_step_"] button {
            background: linear-gradient(135deg, rgba(128, 128, 128, 0.06) 0%, rgba(128, 128, 128, 0.12) 100%) !important;
            border: 2px solid rgba(128, 128, 128, 0.3) !important;
            border-radius: 12px !important;
            min-height: 150px !important;
            padding: 1.25rem 1rem !important;
            box-shadow: 0 2px 8px rgba(0,0,0,0.1) !important;
            transition: all 0.25s ease !important;
            display: flex !important;
            flex-direction: column !important;
            justify-content: center !important;
            align-items: center !important;
            text-align: center !important;
        }
        [class*="st-key-flow_step_"] button p {
            color: var(--fa-accent) !important;
            font-size: 1.15rem !important;
            font-weight: 700 !important;
            margin: 0 !important;
        }
        [class*="st-key-flow_step_"] button::after {
            display: block;
            font-size: 0.9rem;
            font-weight: 500;
            color: #6c757d;
            margin-top: 0.4rem;
        }
        [class*="st-key-flow_step_"] button:hover {
            background: linear-gradient(135deg, #29379b 0%, #1e2a7a 100%) !important;
            border-color: #29379b !important;
            transform: translateY(-3px) !important;
            box-shadow: 0 8px 24px rgba(41, 55, 155, 0.3) !important;
        }
        [class*="st-key-flow_step_"] button:hover p,
        [class*="st-key-flow_step_"] button:hover::after {
            color: white !important;
        }
        .flow-arrow {
            text-align: center;
            font-size: 2rem;
            color: var(--fa-accent);
            opacity: 0.6;
        }
        @media (max-width: 640px) {
            .flow-arrow { transform: rotate(90deg); }
        }
        """
        + "\n".join(css)
        + "</style>",
        unsafe_allow_html=True,
    )

def render_whats_new():
    """Render the release highlights below the workflow buttons."""
    spacer1, center_col, spacer2 = st.columns([1, 6, 1])
    with center_col:
        st.markdown(
            """
            <div class="whats-new">
                <h3>✨ What's new in FLASHApp 2.0 (alpha 10)</h3>
                <div class="whats-new-columns">
                    <div>
                        <h4>Explore results your way</h4>
                        <ul>
                            <li>Switch between ready-made views or make your own; each keeps its zoom and selection.</li>
                            <li>Compare runs in Quantification, Identification and QC views, with the QC of a label-free or TMT/iTRAQ comparison.</li>
                            <li>New: 3D MS1 map, precursor spectrum, isotope envelope and mass shifts of identified proteoforms.</li>
                        </ul>
                    </div>
                    <div>
                        <h4>Sequences and identifications</h4>
                        <ul>
                            <li>Edit an identified proteoform, paste or copy ProForma, and pick variable modifications from a library.</li>
                            <li>Matched fragment ions are marked on the raw and deconvolved spectra.</li>
                            <li>Add identification to a deconvolved run without deconvolving it again, or match a mass list against a sequence.</li>
                        </ul>
                    </div>
                </div>
            </div>
            <style>
            .whats-new {
                background: linear-gradient(135deg, rgba(128, 128, 128, 0.06) 0%, rgba(128, 128, 128, 0.12) 100%);
                border: 1px solid rgba(128, 128, 128, 0.3);
                border-left: 5px solid var(--fa-accent);
                border-radius: 12px;
                padding: 1.5rem 2rem;
                margin: 1.5rem 0 0.5rem 0;
                box-shadow: 0 2px 8px rgba(0,0,0,0.05);
            }
            .whats-new h3 {
                color: var(--fa-accent);
                font-size: 1.4rem;
                font-weight: 700;
                margin: 0 0 1rem 0;
                padding: 0;
            }
            .whats-new-columns {
                display: flex;
                flex-wrap: wrap;
                gap: 1rem 2.5rem;
            }
            .whats-new-columns > div {
                flex: 1 1 260px;
            }
            .whats-new h4 {
                color: var(--fa-accent);
                font-size: 1.05rem;
                font-weight: 600;
                margin: 0 0 0.5rem 0;
                padding: 0;
            }
            .whats-new ul {
                margin: 0;
                padding-left: 1.2rem;
                color: var(--fa-body-text);
            }
            .whats-new li {
                margin-bottom: 0.35rem;
            }
            </style>
            """,
            unsafe_allow_html=True,
        )

# Main execution
def main():
    """Main function to render the quickstart page."""
    # Inject custom CSS
    inject_theme_colors()
    inject_workflow_button_css()
    
    # Render main sections: title and logo, the downloads, the workflow, what's new
    render_workflow_selection()

    render_download_section()

    render_workflow_flow()

    render_whats_new()

    render_contributors()
    

main()
