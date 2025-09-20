"""
Installation Check Script

This script verifies that all required dependencies are installed and functioning correctly.
"""
import sys
import subprocess
import pkg_resources
import importlib
from pathlib import Path

# Required packages and their minimum versions
REQUIRED_PACKAGES = {
    'streamlit': '1.32.0',
    'pandas': '2.1.0',
    'numpy': '1.24.3',
    'plotly': '5.17.0',
    'biopython': '1.81',
    'py3dmol': '2.0.0',
    'stmol': '0.0.6',
    'scikit-learn': '1.3.0',
    'scipy': '1.11.0',
    'matplotlib': '3.7.0',
    'seaborn': '0.12.0',
    'freesasa': '2.1.0',
    'rdkit': '2023.3.0',
    'openbabel': '3.1.0',
    'pytest': '7.4.0'
}

def check_python_version():
    """Check if Python version is compatible."""
    print("\n🔍 Checking Python version...")
    if sys.version_info < (3, 8):
        print(f"❌ Python 3.8 or higher is required. Current version: {sys.version}")
        return False
    print(f"✅ Python {sys.version.split()[0]} is compatible")
    return True

def check_package_installation():
    """Check if all required packages are installed."""
    print("\n🔍 Checking package installations...")
    missing_packages = []
    wrong_version = []
    
    for package, required_version in REQUIRED_PACKAGES.items():
        try:
            installed_version = pkg_resources.get_distribution(package).version
            if pkg_resources.parse_version(installed_version) < pkg_resources.parse_version(required_version):
                wrong_version.append((package, installed_version, required_version))
        except pkg_resources.DistributionNotFound:
            missing_packages.append((package, required_version))
    
    if missing_packages:
        print("❌ Missing packages:")
        for package, version in missing_packages:
            print(f"  - {package} (required: >={version})")
    
    if wrong_version:
        print("⚠️  Outdated packages:")
        for package, installed, required in wrong_version:
            print(f"  - {package}: {installed} (required: >={required})")
    
    if not missing_packages and not wrong_version:
        print("✅ All required packages are installed with compatible versions")
        return True
    
    return False

def check_imports():
    """Check if all required modules can be imported."""
    print("\n🔍 Checking module imports...")
    failed_imports = []
    
    for package in REQUIRED_PACKAGES.keys():
        try:
            importlib.import_module(package.split('<')[0].strip())
            print(f"✅ {package} imports successfully")
        except ImportError as e:
            failed_imports.append((package, str(e)))
    
    if failed_imports:
        print("❌ Failed to import:")
        for package, error in failed_imports:
            print(f"  - {package}: {error}")
        return False
    
    return True

def check_directory_structure():
    """Verify that the required directory structure exists."""
    print("\n🔍 Checking directory structure...")
    required_dirs = [
        'data', 'targets', 'structures', 
        'pockets', 'ligands', 'nanodelivery', 'notebooks'
    ]
    missing_dirs = []
    
    for dir_name in required_dirs:
        if not Path(dir_name).exists():
            missing_dirs.append(dir_name)
    
    if missing_dirs:
        print("❌ Missing directories:")
        for dir_name in missing_dirs:
            print(f"  - {dir_name}")
        return False
    
    print("✅ All required directories exist")
    return True

def check_pdb_utils():
    """Test the PDB utilities."""
    print("\n🔍 Testing PDB utilities...")
    try:
        from pdb_utils import visualize_pdb, extract_ligands, calculate_sasa
        print("✅ PDB utilities import successfully")
        return True
    except Exception as e:
        print(f"❌ Error importing PDB utilities: {str(e)}")
        return False

def run_tests():
    """Run basic tests."""
    print("\n🔍 Running basic tests...")
    try:
        import pytest
        result = subprocess.run([sys.executable, "-m", "pytest", "-v"], 
                              capture_output=True, text=True)
        print(result.stdout)
        if result.returncode != 0:
            print("❌ Tests failed")
            return False
        print("✅ All tests passed")
        return True
    except Exception as e:
        print(f"❌ Error running tests: {str(e)}")
        return False

def main():
    """Main function to run all checks."""
    print("🚀 Starting installation check...")
    
    checks = {
        "Python Version": check_python_version(),
        "Package Installation": check_package_installation(),
        "Module Imports": check_imports(),
        "Directory Structure": check_directory_structure(),
        "PDB Utilities": check_pdb_utils(),
        "Tests": run_tests()
    }
    
    # Summary
    print("\n📊 Check Summary:")
    for check_name, status in checks.items():
        status_icon = "✅" if status else "❌"
        print(f"{status_icon} {check_name}")
    
    if all(checks.values()):
        print("\n🎉 All checks passed! Your environment is ready to use.")
        print("You can now run the application with: streamlit run app.py")
        return 0
    else:
        print("\n❌ Some checks failed. Please address the issues above.")
        return 1

if __name__ == "__main__":
    sys.exit(main())
