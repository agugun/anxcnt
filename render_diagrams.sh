#!/bin/bash
# AXSCNT Diagram Renderer
# Uses the local D2 installation to update all SVGs

D2_EXEC="./d2_bin/bin/d2"

if [ ! -f "$D2_EXEC" ]; then
    echo "Error: D2 executable not found at $D2_EXEC"
    echo "Please ensure you are running this from the repository root."
    exit 1
fi

echo "🎨 Rendering AXSCNT Architecture Diagrams..."

for d2_file in doc/diagrams/*.d2; do
    svg_file="${d2_file%.d2}.svg"
    echo "  -> Rendering $d2_file to $svg_file"
    $D2_EXEC "$d2_file" "$svg_file"
done

echo "✅ All diagrams updated successfully."
