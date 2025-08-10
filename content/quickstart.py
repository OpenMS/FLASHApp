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
            height: 240px;
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
            font-size: 3.5rem;
            margin-bottom: 1rem;
            transition: transform 0.3s ease;
        }
        
        .workflow-title {
            font-size: 1.5rem;
            font-weight: 700;
            color: #29379b;
            margin-bottom: 0.5rem;
            transition: color 0.3s ease;
        }
        
        .workflow-subtitle {
            font-size: 1rem;
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
        
        /* Responsive design */
        @media (max-width: 768px) {
            .workflow-button {
                height: 200px;
                padding: 1.5rem 1rem;
                margin-bottom: 1rem;
            }
            
            .workflow-emoji {
                font-size: 2.5rem;
            }
            
            .workflow-title {
                font-size: 1.25rem;
            }
            
            .workflow-subtitle {
                font-size: 0.9rem;
            }
            
            .hero-title {
                font-size: 2rem;
            }
            
            .hero-subtitle {
                font-size: 1.1rem;
            }
        }
        
        @media (max-width: 480px) {
            .workflow-button {
                height: 180px;
                padding: 1rem;
            }
            
            .workflow-emoji {
                font-size: 2rem;
            }
            
            .workflow-title {
                font-size: 1.1rem;
            }
            
            .hero-title {
                font-size: 1.75rem;
            }
        }
        </style>
        """,
        unsafe_allow_html=True,
    )

def create_workflow_button(emoji, title, subtitle, workflow_key, page_name):
    """Create a workflow selection button with navigation."""
    button_html = f"""
    <div class="workflow-button" onclick="location.href='#{workflow_key}-{page_name}'">
        <div class="workflow-emoji">{emoji}</div>
        <div class="workflow-title">{title}</div>
        <div class="workflow-subtitle">{subtitle}</div>
    </div>
    """
    
    if st.markdown(button_html, unsafe_allow_html=True):
        return True
    return False

def create_navigation_button(emoji, title, subtitle, page_path):
    """Create a workflow button that navigates to the specified page."""
    button_html = f"""
    <div class="workflow-button">
        <div class="workflow-emoji">{emoji}</div>
        <div class="workflow-title">{title}</div>
        <div class="workflow-subtitle">{subtitle}</div>
    </div>
    """
    
    # Display the styled button
    st.markdown(button_html, unsafe_allow_html=True)
    
    # Create a button for navigation
    if st.button(f"Start {title}", key=f"{title.lower()}_nav_btn", use_container_width=True, type="primary"):
        st.switch_page(page_path)

def render_workflow_selection():
    """Render the main workflow selection section."""
    # Hero section
    st.markdown(
        """
        <div class="hero-section">
            <h1 class="hero-title">👋 Quick Start</h1>
            <p class="hero-subtitle">FLASHApp: Choose Your Analysis Workflow</p>
        </div>
        """,
        unsafe_allow_html=True,
    )
    
    # Main workflow selection buttons in 3-column layout
    col1, col2, col3 = st.columns(3)
    
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
    """Render the enhanced Windows download section."""
    if Path("OpenMS-App.zip").exists():
        st.markdown(
            """
            <div class="download-section">
                <h3>🪟 Download for Windows</h3>
                <p>Get the standalone Windows application for offline use</p>
            </div>
            """,
            unsafe_allow_html=True,
        )
        
        # Center the download button
        col1, col2, col3 = st.columns([1, 2, 1])
        with col2:
            with open("OpenMS-App.zip", "rb") as file:
                st.download_button(
                    label="📥 Download for Windows",
                    data=file,
                    file_name="OpenMS-App.zip",
                    mime="archive/zip",
                    type="primary",
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
    
    v_space(2)
    
    render_enhanced_download_section()
    
    render_footer()

# Execute main function
if __name__ == "__main__":
    main()
else:
    main()
