'use strict';

/**
 * HRP Landing Tab Plugin
 * This plugin adds an HRP tab to the index page
 */

(function() {
    // Wait for DOM to be ready
    const initHRPTab = () => {
        const pluginContainer = document.getElementById('hrp_landing_tab_html_c');
        
        if (!pluginContainer) {
            console.error('HRP plugin container not found');
            return;
        }

        // Find the tab button and content in the plugin container
        const tabButton = pluginContainer.querySelector('#hrp-tab-button');
        const tabContent = pluginContainer.querySelector('#hrp-tab-content');
        
        if (!tabButton || !tabContent) {
            console.error('HRP tab elements not found in plugin container');
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

        // Add click handler for the HRP tab button
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

        console.log('HRP landing tab plugin initialized successfully');
    };

    // Initialize when DOM is ready
    if (document.readyState === 'loading') {
        document.addEventListener('DOMContentLoaded', initHRPTab);
    } else {
        // DOM is already ready, initialize immediately
        initHRPTab();
    }
})();
