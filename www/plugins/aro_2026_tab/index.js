'use strict';

/**
 * ARO 2026 Tab Plugin
 * This plugin adds an ARO 2026 workshop tab to the index page
 */

(function() {
    // Wait for DOM to be ready
    const initARO2026Tab = () => {
        const pluginContainer = document.getElementById('aro_2026_tab_html_c');
        
        if (!pluginContainer) {
            console.error('ARO 2026 plugin container not found');
            return;
        }

        // Find the tab button and content in the plugin container
        const tabButton = pluginContainer.querySelector('#aro2026-tab-button');
        const tabContent = pluginContainer.querySelector('#aro2026-tab-content');
        
        if (!tabButton || !tabContent) {
            console.error('ARO 2026 tab elements not found in plugin container');
            return;
        }

        // Find the tabs list and content list in the main page
        const tabsList = document.querySelector('#search-tabs-c .tabs ul.tabs-content');
        const tabsContentList = document.querySelector('#search-tabs-c .tabs-content > ul');
        
        if (!tabsList || !tabsContentList) {
            console.error('Index page tab containers not found');
            return;
        }

        // Move the tab button to the tabs list
        tabsList.appendChild(tabButton);
        
        // Move the tab content to the tabs content list
        tabsContentList.appendChild(tabContent);

        // Add click handler for the ARO 2026 tab button
        const tabLink = tabButton.querySelector('a');
        if (tabLink) {
            tabLink.addEventListener('click', (event) => {
                const tab_id = tabButton.dataset.tabId;

                // Update active tab button
                document.querySelectorAll('#search-tabs-c .tabs li').forEach((element) => {
                    if (element.dataset.tabId === tab_id) {
                        element.classList.add('is-active');
                    } else {
                        element.classList.remove('is-active');
                    }
                });

                // Update active tab content
                document.querySelectorAll('#search-tabs-c .tabs-content li').forEach((element) => {
                    if (element.dataset.tabId === tab_id) {
                        element.classList.add('is-active');
                    } else {
                        element.classList.remove('is-active');
                    }
                });
            });
        }

        // Clean up - remove the plugin container since we've moved its contents
        pluginContainer.remove();
    };

    // Initialize when DOM is ready
    if (document.readyState === 'loading') {
        document.addEventListener('DOMContentLoaded', initARO2026Tab);
    } else {
        initARO2026Tab();
    }
})();
