"""Configuration settings for the Drug Target Platform."""
import os
from pathlib import Path
from typing import Dict, Any, Optional
import yaml

# Base directory
BASE_DIR = Path(__file__).parent.resolve()

# Default configuration
DEFAULT_CONFIG = {
    "data": {
        "input_dir": str(BASE_DIR / "data"),
        "output_dir": str(BASE_DIR / "output"),
        "allowed_extensions": [".csv", ".xlsx", ".xls", ".tsv"],
    },
    "analysis": {
        "default_min_logfc": 1.0,
        "default_max_pval": 0.05,
        "top_n_targets": 50,
    },
    "structure": {
        "af2_dir": str(BASE_DIR / "af2_models"),
        "max_sequence_length": 2500,
        "min_confidence": 70.0,
    },
    "docking": {
        "vina_executable": "vina",
        "prepare_receptor": "prepare_receptor",
        "prepare_ligand": "prepare_ligand",
        "box_size": 20.0,  # Angstroms
        "exhaustiveness": 8,
        "num_modes": 9,
    },
    "nanocarrier": {
        "default_carrier": "liposome",
        "default_ligand": "antibody",
        "default_trigger": "pH",
    },
    "logging": {
        "level": "INFO",
        "file": str(BASE_DIR / "logs" / "app.log"),
        "max_size_mb": 10,
        "backup_count": 3,
    },
}

class Config:
    """Configuration manager for the application."""
    
    _instance = None
    
    def __new__(cls):
        if cls._instance is None:
            cls._instance = super(Config, cls).__new__(cls)
            cls._instance._config = DEFAULT_CONFIG.copy()
            cls._instance._load_config()
        return cls._instance
    
    def _load_config(self):
        """Load configuration from file if it exists."""
        config_path = BASE_DIR / "config.yaml"
        if config_path.exists():
            try:
                with open(config_path, 'r') as f:
                    file_config = yaml.safe_load(f) or {}
                self._deep_update(self._config, file_config)
            except Exception as e:
                print(f"Warning: Could not load config file: {e}")
    
    def _deep_update(self, original: Dict[Any, Any], update: Dict[Any, Any]) -> None:
        """Recursively update a dictionary."""
        for key, value in update.items():
            if key in original and isinstance(original[key], dict) and isinstance(value, dict):
                self._deep_update(original[key], value)
            else:
                original[key] = value
    
    def get(self, key: str, default: Any = None) -> Any:
        """Get a configuration value by dot notation."""
        keys = key.split('.')
        value = self._config
        
        try:
            for k in keys:
                value = value[k]
            return value
        except (KeyError, TypeError):
            return default
    
    def set(self, key: str, value: Any) -> None:
        """Set a configuration value by dot notation."""
        keys = key.split('.')
        config = self._config
        
        for k in keys[:-1]:
            if k not in config or not isinstance(config[k], dict):
                config[k] = {}
            config = config[k]
        
        config[keys[-1]] = value
    
    def save(self) -> bool:
        """Save the current configuration to file."""
        config_path = BASE_DIR / "config.yaml"
        try:
            with open(config_path, 'w') as f:
                yaml.dump(self._config, f, default_flow_style=False)
            return True
        except Exception as e:
            print(f"Error saving config: {e}")
            return False

# Global configuration instance
config = Config()

def get_path(*path_parts: str) -> Path:
    """Get a path relative to the base directory."""
    return Path(BASE_DIR, *path_parts)

def ensure_dir(path: Path) -> Path:
    """Ensure a directory exists and return its path."""
    path.mkdir(parents=True, exist_ok=True)
    return path
