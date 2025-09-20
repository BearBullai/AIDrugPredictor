"""
Install Missing Dependencies

This script installs all required Python packages for the project.
"""
import sys
import subprocess
import importlib

def install_package(package):
    """Install a Python package using pip."""
    try:
        subprocess.check_call([sys.executable, "-m", "pip", "install", package])
        print(f"✅ Successfully installed {package}")
        return True
    except subprocess.CalledProcessError as e:
        print(f"❌ Failed to install {package}: {e}")
        return False

def main():
    """Main function to install missing dependencies."""
    print("🚀 Starting dependency installation...\n")
    
    # List of packages to install
    packages = [
        'biopython>=1.81',
        'py3dmol>=2.0.0',
        'scikit-learn>=1.3.0',
        'freesasa>=2.1.0',
        'openbabel>=3.1.0',
        'pytest>=7.4.0'
    ]
    
    # Install each package
    success = True
    for package in packages:
        print(f"Installing {package}...")
        if not install_package(package):
            success = False
    
    # Special handling for OpenBabel on Windows
    if sys.platform.startswith('win'):
        print("\nℹ️  For Windows users: If OpenBabel installation failed, "
              "please download and install it manually from "
              "https://anaconda.org/conda-forge/openbabel")
    
    # Verify installations
    print("\n🔍 Verifying installations...")
    for package in packages:
        pkg_name = package.split('>=')[0].split('<')[0].strip()
        try:
            importlib.import_module(pkg_name)
            print(f"✅ {pkg_name} is installed and importable")
        except ImportError as e:
            print(f"❌ {pkg_name} import failed: {e}")
            success = False
    
    if success:
        print("\n🎉 All dependencies installed successfully!")
        print("You can now run the application with: streamlit run app.py")
    else:
        print("\n⚠️  Some dependencies failed to install. "
              "Please check the error messages above and install them manually.")
    
    return 0 if success else 1

if __name__ == "__main__":
    sys.exit(main())
