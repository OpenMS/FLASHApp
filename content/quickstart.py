import streamlit as st

from src.common.common import page_setup, themed_asset, v_space
from src.common.installers import installers_available, render_installer_buttons

page_setup(page="main")

def render_download_box(container_key, title, description, icon="", subtitle=""):
    """Render a styled download box container with consistent styling."""
    with st.container(key=container_key):
        st.markdown(
            f"""
            <h4 style="color: #6c757d; margin-bottom: 0.75rem; font-size: 1.1rem; font-weight: 600;">
                {icon} {title}
            </h4>
            <p style="color: #6c757d; margin-bottom: 1rem;">
                {description}
            </p>
            """,
            unsafe_allow_html=True,
        )
        
        # Return context for button placement
        return st.columns([1, 2, 1])

def render_desktop_download_box():
    """Render the desktop app download box: one button per platform."""
    if not installers_available():
        return False

    with st.container(key="desktop_download_container"):
        st.markdown(
            """
            <h4 style="color: #6c757d; margin-bottom: 0.75rem; font-size: 1.1rem; font-weight: 600;">
                FLASHApp for your computer
            </h4>
            <p style="color: #6c757d; margin-bottom: 1rem;">
                FLASHApp is best enjoyed online but you can download an offline version for Windows, macOS (Apple Silicon) and Linux below.
            </p>
            """,
            unsafe_allow_html=True,
        )
        render_installer_buttons(app_name="FLASHApp")
        st.markdown(
            """
            <div style="text-align: center; margin-top: 1rem; color: #6c757d;">
                No admin rights are needed on Windows. On macOS, drag FLASHApp from the
                disk image to Applications. On Linux, make the AppImage executable and start it.<br>
                Even offline, it's still a web app - just packaged so you can use it without an internet connection.
            </div>
            """,
            unsafe_allow_html=True,
        )

    return True

def render_flash_cli_tools_box():
    """Render the FLASH* command line tools box."""
    cols = render_download_box(
        container_key="flash_cli_tools_container",
        title="FLASH* Command Line Tools",
        description="Access FLASH* tools (FLASHDeconv, FLASHTnT, etc.) through the command line interface as part of OpenMS.",
        icon=""
    )
    
    # Center the link button using Streamlit button with custom styling
    with cols[1]:

        st.link_button(
            "📥 Download Command Line Tools",
            "https://abibuilder.cs.uni-tuebingen.de/archive/openms/OpenMSInstaller/experimental/FVdeploy/",
            use_container_width=True,
            type="secondary",   # matches the solid/primary button look
        )

    
    # Add subtitle text
    st.markdown(
        """
        <div style="text-align: center; margin-top: 1rem; color: #6c757d;">
            FLASH* tools are a part of the OpenMS TOPP tools and can be used via command line interface for advanced custom automated workflows and scripting.
        </div>
        """,
        unsafe_allow_html=True,
    )

def apply_download_box_styling(container_key):
    """Apply consistent styling to a download box container."""
    st.markdown(
        f"""
        <style>
        .st-key-{container_key} {{
            background: linear-gradient(135deg, rgba(128, 128, 128, 0.06) 0%, rgba(128, 128, 128, 0.12) 100%) !important;
            border: 1px solid rgba(128, 128, 128, 0.3) !important;
            border-radius: 8px !important;
            padding: 1.5rem !important;
            margin: 1rem 0 !important;
            text-align: center !important;
            box-shadow: 0 2px 4px rgba(0,0,0,0.05) !important;
        }}
        
        .st-key-{container_key} > div {{
            background: transparent !important;
        }}
        </style>
        """,
        unsafe_allow_html=True,
    )

def render_download_boxes_responsive():
    """Render download boxes in a responsive layout."""
    # Collect available boxes
    available_boxes = []
    
    # Check if a desktop installer is available
    if installers_available():
        available_boxes.append("desktop")
    
    # Always include FLASH* CLI tools
    available_boxes.append("flash_cli")
    
    num_boxes = len(available_boxes)
    
    if num_boxes == 0:
        return
    
    # Section header
    st.markdown(
        """
        <h4 style="color: #6c757d; margin-bottom: 1.5rem; font-size: 1.3rem; font-weight: 600; text-align: center;">
            Want to use FLASHApp offline?
        </h4>
        """,
        unsafe_allow_html=True,
    )

    spacer1, center_col, spacer2 = st.columns([1, 3, 1])
    # Render boxes based on count
    if num_boxes == 1:
        # Single box: centered
        with center_col:
            render_flash_cli_tools_box()
            apply_download_box_styling("flash_cli_tools_container")
                
    elif num_boxes == 2:
        with center_col:
            render_desktop_download_box()
            apply_download_box_styling("desktop_download_container")
            render_flash_cli_tools_box()
            apply_download_box_styling("flash_cli_tools_container")


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

def create_navigation_button(emoji, title, subtitle, page_path):
    """Create a functional workflow button that navigates to the specified page."""
    
    # Create unique key for this button
    button_key = f"{title.lower().replace(' ', '_')}_nav_btn"
    
    # Create the button with custom styling applied via CSS classes
    button_label = f"{emoji} {title}"
    
    # Use Streamlit's button with custom styling
    if st.button(
        label=button_label,
        key=button_key,
        help=f"Navigate to {title} - {subtitle}",
        use_container_width=True,
        type="primary"
    ):
        st.switch_page(page_path)
    
    # Apply custom CSS styling using the key-based selector approach
    st.markdown(
        f"""
        <style>
        /* Target the specific button using key-based selector */
        .st-key-{button_key} button {{
            background: linear-gradient(135deg, rgba(128, 128, 128, 0.06) 0%, rgba(128, 128, 128, 0.12) 100%) !important;
            border: 2px solid rgba(128, 128, 128, 0.3) !important;
            border-radius: 12px !important;
            padding: 2rem 1.5rem !important;
            height: max(280px, 20vh) !important;
            min-height: max(280px, 20vh) !important;
            box-shadow: 0 2px 8px rgba(0,0,0,0.1) !important;
            transition: all 0.3s ease !important;
            margin-bottom: 1rem !important;
            color: var(--fa-accent) !important;
            font-size: clamp(1.25rem, 2.5vw, 1.5rem) !important;
            white-space: nowrap !important;
            overflow: hidden !important;
            text-overflow: ellipsis !important;
            font-weight: 700 !important;
            display: flex !important;
            flex-direction: column !important;
            justify-content: center !important;
            align-items: center !important;
            text-align: center !important;
            line-height: 1.4 !important;
            width: 100% !important;
        }}
        
        .st-key-{button_key} button p {{
            color: var(--fa-accent) !important;
            font-size: clamp(1.25rem, 2.5vw, 1.5rem) !important;
            white-space: nowrap !important;
            overflow: hidden !important;
            text-overflow: ellipsis !important;
            font-weight: 700 !important;
            margin: 0 !important;
        }}
        
        .st-key-{button_key} button:hover {{
            background: linear-gradient(135deg, #29379b 0%, #1e2a7a 100%) !important;
            border-color: #29379b !important;
            transform: translateY(-4px) !important;
            box-shadow: 0 8px 24px rgba(41, 55, 155, 0.3) !important;
        }}
        
        .st-key-{button_key} button:hover p {{
            color: white !important;
        }}
        
        .st-key-{button_key} button:active {{
            background: linear-gradient(135deg, #1e2a7a 0%, #162159 100%) !important;
            border-color: #1e2a7a !important;
            transform: translateY(-2px) scale(0.98) !important;
            box-shadow: 0 4px 12px rgba(41, 55, 155, 0.4) !important;
        }}
        
        .st-key-{button_key} button:active p {{
            color: white !important;
        }}
        
        .st-key-{button_key} button:focus {{
            background: linear-gradient(135deg, #29379b 0%, #1e2a7a 100%) !important;
            border-color: #29379b !important;
        }}
        
        .st-key-{button_key} button:focus p {{
            color: white !important;
        }}
        
        /* Add subtitle styling using pseudo-element */
        .st-key-{button_key} button::after {{
            content: "{subtitle}";
            display: block;
            font-size: clamp(0.875rem, 1.5vw, 1rem) !important;
            white-space: nowrap !important;
            overflow: hidden !important;
            text-overflow: ellipsis !important;
            font-weight: 500 !important;
            color: #6c757d !important;
            margin-top: 0.5rem !important;
        }}
        
        .st-key-{button_key} button:hover::after {{
            color: #e9ecef !important;
        }}
        
        .st-key-{button_key} button:active::after {{
            color: #e9ecef !important;
        }}
        
        .st-key-{button_key} button:focus::after {{
            color: #e9ecef !important;
        }}
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
    
    # Main workflow selection buttons with centered, compact layout
    # Use spacing columns to center buttons and prevent wide-screen spreading
    spacer1, col1, col2, col3, spacer2 = st.columns([1, 2, 2, 2, 1], gap="small")
    
    with col1:
        create_navigation_button(
            "⚡️",
            "Deconvolution",
            "FLASHDeconv",
            "content/FLASHDeconv/FLASHDeconvWorkflow.py"
        )
    
    with col2:
        create_navigation_button(
            "🧨",
            "Identification",
            "FLASHTnT",
            "content/FLASHTnT/FLASHTnTWorkflow.py"
        )
    
    with col3:
        create_navigation_button(
            "📊",
            "Quantification",
            "FLASHQuant",
            "content/FLASHQuant/FLASHQuantFileUpload.py"
        )

def render_whats_new():
    """Render the release highlights below the workflow buttons."""
    spacer1, center_col, spacer2 = st.columns([1, 6, 1])
    with center_col:
        st.markdown(
            """
            <div class="whats-new">
                <h3>✨ What's new in FLASHApp 1.1.0</h3>
                <div class="whats-new-columns">
                    <div>
                        <h4>Faster and smoother</h4>
                        <ul>
                            <li>Results open much faster, even for large datasets.</li>
                            <li>Your analysis finishes sooner after the tools have run.</li>
                            <li>The viewer stays responsive while you click through scans and proteins.</li>
                        </ul>
                    </div>
                    <div>
                        <h4>Better identification</h4>
                        <ul>
                            <li>Updated to the latest FLASHDeconv and FLASHTnT.</li>
                            <li>More fragment ions are identified, matched directly on the raw spectra.</li>
                            <li>New one-click settings for Lysate and Native samples.</li>
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

def render_enhanced_download_section():
    """Render the enhanced download section with modular components."""
    # The offline downloads are pointless inside the offline app itself.
    if st.session_state.get("location") == "local":
        return
    # Add spacing
    st.markdown("<br><hr><br>", unsafe_allow_html=True)
    
    # Render download boxes using responsive layout
    render_download_boxes_responsive()

# Main execution
def main():
    """Main function to render the quickstart page."""
    # Inject custom CSS
    inject_theme_colors()
    inject_workflow_button_css()
    
    # Render main sections
    render_workflow_selection()

    render_whats_new()

    render_enhanced_download_section()
    

main()
