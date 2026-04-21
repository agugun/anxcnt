import os
import importlib
import json

class ModuleRegistry:
    def __init__(self, modules_dir="web/modules"):
        self.modules_dir = modules_dir
        self.modules = {}

    def discover_modules(self):
        """Scans the modules directory for subfolders with layout.py"""
        if not os.path.exists(self.modules_dir):
            return {}

        for entry in os.listdir(self.modules_dir):
            module_path = os.path.join(self.modules_dir, entry)
            if os.path.isdir(module_path):
                # Look for layout.py which is the entry point for custom UI
                if "layout.py" in os.listdir(module_path):
                    module_name = entry
                    try:
                        # Dynamic import: web.modules.reservoir.layout
                        pkg_path = f"web.modules.{module_name}.layout"
                        mod = importlib.import_module(pkg_path)
                        
                        # Use manifest.json if present, else use folder name
                        manifest_path = os.path.join(module_path, "manifest.json")
                        title = module_name.capitalize()
                        icon = "bi-box"
                        description = "Domain-specific analyzer"
                        if os.path.exists(manifest_path):
                            with open(manifest_path, 'r') as f:
                                manifest = json.load(f)
                                title = manifest.get("title", title)
                                icon = manifest.get("icon", icon)
                                description = manifest.get("description", description)
                        
                        self.modules[module_name] = {
                            "id": module_name,
                            "title": title,
                            "icon": icon,
                            "description": description,
                            "layout": getattr(mod, "layout", None),
                            "path": pkg_path
                        }
                    except Exception as e:
                        print(f"Error loading module {module_name}: {e}")
        
        return self.modules

    def get_module_list(self):
        return [
            {
                "label": info["title"], 
                "value": mod_id, 
                "icon": info.get("icon", "bi-box"),
                "description": info.get("description", "")
            } 
            for mod_id, info in self.modules.items()
        ]
