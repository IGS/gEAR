
// Class for the faceted search widget
// SAdkins - I let Github Copilot write this class for me (with minor adjustments)
class FacetWidget {
    constructor({aggregations, filters, onFilterChange, facetContainer, selectedFacetsTags}) {
        this.aggregations = aggregations || {}; // Counts of all categories (filter)
        this.filters = filters || {};   // Selected categories and values
        this.onFilterChange = onFilterChange || {}; // Callback function when a filter is changed
        this.facetContainer =  facetContainer || document.getElementById('facet_c');
        this.selectedFacetsTags = selectedFacetsTags || document.getElementById('selected_facets_tags');
        this.init();
    }

    init() {
        // Remove any existing facet widget
        const facetWidget = document.getElementById('facet-widget');
        if (facetWidget !== null) {
            facetWidget.remove();
        }
        if (this.aggregations.length > 0 ) {
            this.createFacetWidget();
        } else {
            console.warn('No aggregations to display');
        }
    }

    createFacetWidget() {
        const facetWidget = document.createElement('div');
        facetWidget.id = 'facet-widget';
        facetWidget.className = 'facet-widget content is-small columns';

        this.facetContainer.append(facetWidget);
        this.render();
    }

    updateAggregations(aggregations) {
        this.aggregations = aggregations;
        // Update counts
        for (const filter of this.aggregations) {
            const escapedFilterName = CSS.escape(filter.name);

            for (const item of filter.items) {
                const escapedItemName = CSS.escape(item.name);
                // NOTE: For some reason, CSS.escape() does not escape properly for document.getElementById()
                const filterItemCount = document.querySelector(`#filter-item-${escapedFilterName}-${escapedItemName} .filter-item-count`);
                const filterItem = document.querySelector(`#filter-item-${escapedFilterName}-${escapedItemName}`);
                // show item
                filterItemCount.classList.remove('is-hidden');
                filterItem.classList.remove('is-hidden');
                filterItemCount.innerHTML = item.count;
                // if count is 0, hide item count
                if (item.count === 0) {
                    filterItemCount.classList.add('is-hidden');
                    // if item category is not in the filter category, hide item
                    if (this.filters.hasOwnProperty(filter.name)) {
                        // If item is in the filter category, show 0 item count
                        if (this.filters[filter.name].includes(item.name)) {
                            filterItemCount.classList.remove('is-hidden');
                        }
                    } else {
                        filterItem.classList.add('is-hidden');
                    }

                }
            }
        }
    }

    render() {
        const facetWidget = document.getElementById('facet-widget');
        facetWidget.innerHTML = '';
        const facetWidgetContent = document.createElement('div');
        facetWidgetContent.id = 'facet-widget-content';
        facetWidgetContent.className = 'facet-widget-content column is-two-thirds';
        facetWidget.appendChild(facetWidgetContent);
        this.renderFilters();
    }

    renderFilters() {
        const facetWidgetContent = document.getElementById('facet-widget-content');
        const filters = this.aggregations;
        for (const filter of filters) {
            const filterElement = this.createFilterElement(filter);
            facetWidgetContent.appendChild(filterElement);
        }
        // add event listeners for filter toggle
        const filterHeaders = document.querySelectorAll('.filter-header');
        for (const filterHeader of filterHeaders) {
            filterHeader.addEventListener('click', (event) => {
                this.onFilterHeaderClick(event);
            });
        }

    }

    createFilterElement(filter) {
        const filterElement = document.createElement('div');
        filterElement.id = `filter-${filter.name}`;
        filterElement.className = 'filter';
        filterElement.innerHTML = `
            <div class="filter-header panel-heading is-clickable has-background-primary-dark has-text-primary-light">
                <span class="filter-name has-text-weight-semibold is-underlined">${filter.name}</span>
                <span class="loader is-hidden is-inline-flex"></span>
                <span class="filter-toggle-icon icon is-pulled-right">
                    <i class="mdi mdi-chevron-down"></i>
                </span>
            </div>
        `;
        const filterContent = this.createFilterContent(filter);
        filterElement.appendChild(filterContent);
        return filterElement;
    }

    createFilterContent(filter) {
        const filterContent = document.createElement('div');
        filterContent.id = `filter-content-${filter.name}`;
        filterContent.className = 'filter-content is-hidden has-background-white';
        const filterItems = this.createFilterItems(filter);
        filterContent.appendChild(filterItems);
        return filterContent;
    }

    createFilterItems(filter) {
        const filterItems = document.createElement('div');
        filterItems.id = `filter-items-${filter.name}`;
        filterItems.className = 'filter-items';
        for (const item of filter.items) {
            const filterItem = this.createFilterItem(item, filter.name);
            filterItems.appendChild(filterItem);
        }
        return filterItems;
    }

    createFilterItem(item, filterName) {
        const filterItem = document.createElement('div');
        filterItem.id = `filter-item-${filterName}-${item.name}`;   // some item names can be duplicated across different filters
        filterItem.className = 'filter-item';
        filterItem.innerHTML = `
            <label class="panel-block checkbox">
                <input type="checkbox">
                <span class="filter-item-name has-text-weight-medium">${item.name}</span>
                <span class="filter-item-count tag is-rounded">${item.count}</span>
                <span class="loader is-hidden"></span>
            </label>
        `;

        // Since clicking anything in label will trigger the checkbox, we need to add the event listener to the input
        const inputElt = filterItem.querySelector('input');

        // Check if the filter is already selected
        if (this.filters.hasOwnProperty(filterName) && this.filters[filterName].includes(item.name)) {
            inputElt.checked = true;
        }

        inputElt.addEventListener('click', (event) => {
            this.onFilterItemClick(item, filterName)
        });
        return filterItem;
    }

    onFilterHeaderClick(event) {
        const filterHeader = event.currentTarget;
        const filterContent = filterHeader.nextElementSibling;
        const filterToggleIcon = filterHeader.querySelector('.filter-toggle-icon');
        filterContent.classList.toggle('is-hidden');

        // switch toggle icon to up or down
        if (filterToggleIcon.innerHTML.trim() === '<i class="mdi mdi-chevron-down"></i>') {
            filterToggleIcon.innerHTML = '<i class="mdi mdi-chevron-up"></i>';
        } else {
            filterToggleIcon.innerHTML = '<i class="mdi mdi-chevron-down"></i>';
        }
    }

    async onFilterItemClick(item, filterName) {
        this.toggleFilterItem(item, filterName);
        const escapedFilterName = CSS.escape(filterName);
        const loadingSelector = document.querySelector(`#filter-${escapedFilterName} .loader`);
        loadingSelector.classList.remove("is-hidden")
        await this.onFilterChange(this.filters);
        this.toggleSelectedTag(item.name, filterName);
        loadingSelector.classList.add("is-hidden")

    }

    toggleFilterItem(item, filterName) {
        if (this.filters.hasOwnProperty(filterName)) {
            if (this.filters[filterName].includes(item.name)) {
                this.filters[filterName] = this.filters[filterName].filter((value) => value !== item.name);
            } else {
                this.filters[filterName].push(item.name);
            }
        } else {
            this.filters[filterName] = [item.name];
        }

        // If all items are selected, remove the filter
        if (this.filters[filterName].length === this.aggregations.find((value) => value.name === filterName).items.length) {
            delete this.filters[filterName];
        }

        // If no filters are selected, remove the filter
        if (this.filters[filterName].length === 0) {
            delete this.filters[filterName];
        }
    }

    // Add a tag to the selected filter items list
    toggleSelectedTag(itemName, filterName) {
        // Append the tag if it doesn't exist, otherwise remove it
        const tagId = `selected-tag-${filterName}-${itemName}`
        const selectedTag = document.getElementById(tagId);
        if (selectedTag === null) {
            const tag = document.createElement('span');
            tag.id = tagId;
            tag.className = 'tag is-light is-primary';
            tag.textContent = `${filterName}:${itemName}`;
            this.selectedFacetsTags.appendChild(tag);
            return;
        }
        selectedTag.remove();

    }

}