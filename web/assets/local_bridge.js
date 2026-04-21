/**
 * Local File Bridge for Axscnt Edge-Client
 * Implements the "Philosophy" of Cloud-App + Local-Data
 */

window.axscntBridge = {
    folderHandle: null,

    /**
     * Triggered by the "Connect Local Workspace" button
     */
    async connectWorkspace() {
        try {
            // Native File System Access API
            this.folderHandle = await window.showDirectoryPicker();
            console.log("✅ Local workspace connected:", this.folderHandle.name);
            
            // Notify Dash via a custom signal if needed
            alert(`Connected to: ${this.folderHandle.name}\nLogic is now processing your local disk handle.`);
            
            // In a full implementation, we would register this handle 
            // to a Service Worker to intercept 'fetch' calls to simulation data.
        } catch (err) {
            console.error("❌ Failed to connect workspace:", err);
        }
    },

    /**
     * Audit local folder WITHOUT server traffic
     */
    async auditWorkspace() {
        if (!this.folderHandle) return ["❌ Error: No workspace connected. Please click 'Connect' first."];
        
        let logs = [];
        logs.push(`🔍 Auditing local handle: ${this.folderHandle.name}...`);
        
        try {
            for await (const entry of this.folderHandle.values()) {
                const kind = entry.kind === 'file' ? '📄' : '📁';
                logs.push(`${kind} Found: ${entry.name}`);
            }
            logs.push("✅ Audit Complete. All files scanned locally in-browser.");
            logs.push("🚀 ZERO binary data was sent to server during this scan.");
        } catch (err) {
            logs.push(`❌ Audit Error: ${err.message}`);
        }
        return logs;
    }
};

// Hook into Dash after it loads
document.addEventListener('DOMContentLoaded', () => {
    let checkBtn = setInterval(() => {
        const btn = document.getElementById('btn-connect-local');
        if (btn) {
            btn.onclick = () => window.axscntBridge.connectWorkspace();
            clearInterval(checkBtn);
        }
    }, 500);
});
