from pathlib import Path
import streamlit as st

from src.common.common import page_setup, v_space

page_setup(page="main")

def inject_workflow_button_css():
    """Inject CSS for custom workflow button styling with responsive design."""
    st.markdown(
        """
        <style>
        /* Main workflow button styling */
        .workflow-button {
            background: linear-gradient(135deg, #f8f9fa 0%, #e9ecef 100%);
            border: 2px solid #dee2e6;
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
            color: #29379b;
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
            background: linear-gradient(135deg, #f1f3f4 0%, #e8eaed 100%);
            border: 2px solid #dadce0;
            border-radius: 12px;
            padding: 2rem;
            margin: 2rem 0;
            text-align: center;
        }
        
        .download-section h3 {
            color: #29379b;
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
            color: #29379b;
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
            border-top: 1px solid #dee2e6;
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
            background: linear-gradient(135deg, #f8f9fa 0%, #e9ecef 100%) !important;
            border: 2px solid #dee2e6 !important;
            border-radius: 12px !important;
            padding: 2rem 1.5rem !important;
            height: max(280px, 20vh) !important;
            min-height: max(280px, 20vh) !important;
            box-shadow: 0 2px 8px rgba(0,0,0,0.1) !important;
            transition: all 0.3s ease !important;
            margin-bottom: 1rem !important;
            color: #29379b !important;
            font-size: clamp(1.25rem, 2.5vw, 2rem) !important;
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
            color: #29379b !important;
            font-size: clamp(1.25rem, 2.5vw, 2rem) !important;
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
            font-size: clamp(0.875rem, 1.5vw, 1.25rem) !important;
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
    # Hero section with title on left and OpenMS logo on right
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
        st.image("assets/OpenMS.png", width=200)
    
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

def render_enhanced_download_section():
    if Path("OpenMS-App.zip").exists():
        # Add spacing
        st.markdown("<br><hr><br>", unsafe_allow_html=True)
        
        st.markdown(
            """
            <div style="
                background: linear-gradient(135deg, #f8f9fa 0%, #f1f3f4 100%);
                border: 1px solid #e0e0e0;
                border-radius: 8px;
                padding: 1.5rem;
                margin: 1rem 0;
                text-align: center;
            ">
                <h4 style="color: #6c757d; margin-bottom: 0.75rem; font-size: 1.1rem; font-weight: 600;">
                    Want to use FLASHApp offline?
                </h4>
                <p style="color: #868e96; font-size: 0.9rem; margin-bottom: 1rem;">
                    FLASHApp is best enjoyed online but you can download an offline version for Windows systems below.
                </p>
            </div>
            """,
            unsafe_allow_html=True,
        )
        
        # Center the download button with smaller size
        col1, col2, col3 = st.columns([2, 2, 2])
        with col2:
            with open("OpenMS-App.zip", "rb") as file:
                st.download_button(
                    label="📥 Download ZIP",
                    data=file,
                    file_name="OpenMS-App.zip",
                    mime="archive/zip",
                    type="secondary",
                    use_container_width=True,
                )
        
        st.markdown(
            """
            <div style="text-align: center; margin-top: 1rem; color: #6c757d;">
                Extract the zip file and run the installer (.msi) to install the app.<br>
                Launch using the desktop icon after installation.
            </div>
            """,
            unsafe_allow_html=True,
        )

def render_footer():
    """Render the footer section with new features and OpenMS logo."""
    st.markdown(
        """
        <div class="footer-section">
        """,
        unsafe_allow_html=True,
    )
    
    col1, col2 = st.columns([2, 1])
    
    with col1:
        st.markdown(
            """
            ## ⭐ What's New
            
            **🔄 FLASHViewer is now FLASHApp**
            - Enhanced workflow selection interface
            - Improved navigation and user experience
            - Modern, responsive design for all devices
            
            **🔗 Share & Collaborate**
            - Bookmark your progress with shareable URLs
            - Team collaboration made simple
            - Resume work from any device
            
            **⚡ Performance Improvements**
            - Faster data processing and visualization
            - Optimized memory usage
            - Enhanced stability and reliability
            """
        )
    
    with col2:
        st.image("assets/OpenMS.png", width=300)
    
    st.markdown("</div>", unsafe_allow_html=True)

# Main execution
def main():
    """Main function to render the quickstart page."""
    # Inject custom CSS
    inject_workflow_button_css()
    
    # Render main sections
    render_workflow_selection()
    
    
    render_enhanced_download_section()
    
    render_footer()

main()
