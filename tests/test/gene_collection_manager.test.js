#!/usr/bin/env mocha

/*
    Unit tests for gene_collection_manager.js
*/

// Some setup help by https://publishing-project.rivendellweb.net/testing-front-end-with-mocha-and-playwright-2/

import "mocha" // to set global variables like describe, it, etc.
import { expect } from "playwright/test";

import { browsers, gearBase, login, loginFailure, setupBrowserContext, teardownBrowserContext } from "./helpers.js";

let page;

const gearUrl = `${gearBase}/gene_collection_manager.html`;

/**
 * Mocks the search gene collections API endpoint for testing purposes.
 * @param {Object} page - The page object.
 * @returns {Promise<void>}
 */
const mockSearchGeneCollections = async (page) => {
    await page.route(`${gearBase}/cgi/search_gene_carts.cgi`, async route => {
        const json = {
            "success": 1,
            "problem": "",
            "gene_carts": [
                "{\"id\": 334, \"user_id\": 1, \"gctype\": \"unweighted-list\", \"label\": \"Unweighted test\", \"organism_id\": 1, \"ldesc\": \"This is the unweighted description\", \"share_id\": \"06c1cdd3\", \"is_public\": 1, \"is_domain\": null, \"date_added\": \"2023-11-01 17:10:36\", \"genes\": [\"0610040J01Rik\", \"Zkscan1\", \"1500011K16Rik\", \"1700025G04Rik\", \"Zfp60\", \"1810037I17Rik\", \"2010111I01Rik\", \"2300009A05Rik\", \"Zfp24\", \"Zdhhc2\", \"4931406P16Rik\", \"Ypel2\", \"5031439G07Rik\", \"Wwp1\", \"Wsb2\", \"Wls\", \"Wasf1\", \"Vwc2\", \"Vwa1\", \"Aatk\", \"Vps36\", \"Vps37b\", \"Abca8a\", \"Vldlr\", \"Abhd4\", \"Vcl\", \"Utrn\", \"Acsbg1\", \"Usp1\", \"Adam10\", \"Adam17\", \"Adam23\", \"Adamts17\", \"Adamts5\", \"Ung\", \"Unc119\", \"Adgrg6\", \"Ucp2\", \"Ado\", \"Ado\", \"Uchl1\", \"Adra2c\", \"Afap1l2\", \"Ubl3\", \"Ube2e2\", \"Ube2c\", \"Ahr\", \"Ak3\", \"Ak6\", \"Tyms\", \"Akr1b8\", \"Twf1\", \"Alcam\", \"Alad\", \"Tubg1\", \"Ttyh1\", \"Ttyh2\", \"Tst\", \"Tspo\", \"Angpt2\", \"Ank3\", \"Tspan18\", \"Tspan13\", \"Tspan15\", \"Tspan17\", \"Tsc22d4\", \"Anxa2\", \"Trim2\", \"Trim13\", \"Arfip1\", \"Trappc2\", \"Arhgap19\", \"Arhgap24\", \"Arhgap32\", \"Arhgap39\", \"Arhgef10\", \"Arhgef28\", \"Tpst1\", \"Tppp3\", \"Tpp1\", \"Arid5b\", \"Tpd52\", \"Tox\", \"Arl8a\", \"Arnt2\", \"Top2a\", \"Arrdc3\", \"Art3\", \"Arvcf\", \"Arvcf\", \"Arxes2\", \"Tns3\", \"Asb8\", \"Aspa\", \"Asrgl1\", \"Atad2\", \"Asxl3\", \"Tnik\", \"Atad2\", \"Tnfrsf21\", \"Atg3\", \"Tnfaip6\", \"Tmpo\", \"Atp1b1\", \"Tmod2\", \"Tmem9b\", \"Tmem63b\", \"Atp8a1\", \"Aurkb\", \"Tmem245\", \"B4galt6\", \"Bach1\", \"Bambi\", \"Tmem117\", \"Tm7sf3\", \"Bcas1\", \"Tk1\", \"Tjp1\", \"Bcl7a\", \"Thtpa\", \"Bicd1\", \"Bin3\", \"Birc5\", \"Tgfbr3\", \"Borcs5\", \"Bpnt1\", \"Tex30\", \"Tdrkh\", \"Tead1\", \"Bzw2\", \"Tcp11l1\", \"Tcf19\", \"Tbx2\", \"Cab39l\", \"Cadm1\", \"Cadm4\", \"Camk2d\", \"Tbc1d10a\", \"Tax1bp3\", \"Capg\", \"Tanc2\", \"Capn5\", \"Taldo1\", \"Tagln2\", \"Taf13\", \"Tacc1\", \"Syt11\", \"Syngr1\", \"Cbx6\", \"Svip\", \"Ccdc13\", \"Ccdc28b\", \"Stx7\", \"Strbp\", \"Ccna2\", \"Ccnd1\", \"Ccnb2\", \"Stmn2\", \"Ccser2\"], \"folder_id\": null, \"folder_parent_id\": null, \"folder_label\": null, \"user_name\": \"Test Armstrong\", \"gene_count\": 159, \"organism\": \"Mus musculus\", \"is_owner\": true}",
                "{\"id\": 332, \"user_id\": 1, \"gctype\": \"weighted-list\", \"label\": \"Weighted test\", \"organism_id\": 1, \"ldesc\": \"This is the weighted description\", \"share_id\": \"b5102349\", \"is_public\": 0, \"is_domain\": null, \"date_added\": \"2023-10-24 14:32:04\", \"genes\": [\"0610025J13Rik\", \"0610030E20Rik\", \"0610030E20Rik\", \"0610033M10Rik\", \"0610040B10Rik\", \"1110002L01Rik\", \"1110008L16Rik\", \"1110008P14Rik\", \"1110017D15Rik\", \"1110020A21Rik\", \"1110032F04Rik\"], \"folder_id\": null, \"folder_parent_id\": null, \"folder_label\": null, \"user_name\": \"Test Armstrong\", \"gene_count\": 11, \"organism\": \"Mus musculus\", \"is_owner\": true}",
                "{\"id\": 202, \"user_id\": 2, \"gctype\": \"unweighted-list\", \"label\": \"demo gene collection\", \"organism_id\": 1, \"ldesc\": null, \"share_id\": \"2fe67c37\", \"is_public\": 1, \"is_domain\": null, \"date_added\": \"2022-02-16 14:27:24\", \"genes\": [\"Pde1c\", \"Nav1\", \"Zfp37\", \"Chl1\", \"Mak16\", \"Gab1\", \"Stard9\", \"Farp2\", \"Galnt6\", \"Hjurp\", \"Sept9\", \"Gm14150\", \"Gm16139\", \"Gm35013\"], \"folder_id\": null, \"folder_parent_id\": null, \"folder_label\": null, \"user_name\": \"Test Armstrong\", \"gene_count\": 14, \"organism\": \"Mus musculus\", \"is_owner\": false}",
                //"{\"id\": 28, \"user_id\": 2, \"gctype\": \"unweighted-list\", \"label\": \"test1\", \"organism_id\": 1, \"ldesc\": null, \"share_id\": \"1c9b8f1c\", \"is_public\": 0, \"is_domain\": null, \"date_added\": \"2021-11-23 18:52:36\", \"genes\": [\"arrdc2\", \"kdelr2b\", \"magt1\", \"fosb\", \"aplp2\", \"ak2\", \"szrd1\", \"tmem14ca\", \"cav2\", \"nap1l1\"], \"folder_id\": null, \"folder_parent_id\": null, \"folder_label\": null, \"user_name\": \"Test Armstrong\", \"gene_count\": 10, \"organism\": \"Mus musculus\", \"is_owner\": false}"
            ],
            "pagination": {
                "total_results": 3,
                "current_page": "1",
                "limit": "20",
                "total_pages": 1,
                "next_page": null,
                "prev_page": null
            }
        }
        await route.fulfill({ json });
    });
}

/**
 * Mocks the search gene collections function for a user who is not logged in.
 * @param {Page} page - The page object for the test.
 * @returns {Promise<void>}
 */
const mockSearchGeneCollectionsNotLoggedIn = async (page) => {
    await page.route(`${gearBase}/cgi/search_gene_carts.cgi`, async route => {
        const json = {
            "success": 1,
            "problem": "",
            "gene_carts": [
                "{\"id\": 334, \"user_id\": 1, \"gctype\": \"unweighted-list\", \"label\": \"Unweighted test\", \"organism_id\": 1, \"ldesc\": \"This is the unweighted description\", \"share_id\": \"06c1cdd3\", \"is_public\": 1, \"is_domain\": null, \"date_added\": \"2023-11-01 17:10:36\", \"genes\": [\"0610040J01Rik\", \"Zkscan1\", \"1500011K16Rik\", \"1700025G04Rik\", \"Zfp60\", \"1810037I17Rik\", \"2010111I01Rik\", \"2300009A05Rik\", \"Zfp24\", \"Zdhhc2\", \"4931406P16Rik\", \"Ypel2\", \"5031439G07Rik\", \"Wwp1\", \"Wsb2\", \"Wls\", \"Wasf1\", \"Vwc2\", \"Vwa1\", \"Aatk\", \"Vps36\", \"Vps37b\", \"Abca8a\", \"Vldlr\", \"Abhd4\", \"Vcl\", \"Utrn\", \"Acsbg1\", \"Usp1\", \"Adam10\", \"Adam17\", \"Adam23\", \"Adamts17\", \"Adamts5\", \"Ung\", \"Unc119\", \"Adgrg6\", \"Ucp2\", \"Ado\", \"Ado\", \"Uchl1\", \"Adra2c\", \"Afap1l2\", \"Ubl3\", \"Ube2e2\", \"Ube2c\", \"Ahr\", \"Ak3\", \"Ak6\", \"Tyms\", \"Akr1b8\", \"Twf1\", \"Alcam\", \"Alad\", \"Tubg1\", \"Ttyh1\", \"Ttyh2\", \"Tst\", \"Tspo\", \"Angpt2\", \"Ank3\", \"Tspan18\", \"Tspan13\", \"Tspan15\", \"Tspan17\", \"Tsc22d4\", \"Anxa2\", \"Trim2\", \"Trim13\", \"Arfip1\", \"Trappc2\", \"Arhgap19\", \"Arhgap24\", \"Arhgap32\", \"Arhgap39\", \"Arhgef10\", \"Arhgef28\", \"Tpst1\", \"Tppp3\", \"Tpp1\", \"Arid5b\", \"Tpd52\", \"Tox\", \"Arl8a\", \"Arnt2\", \"Top2a\", \"Arrdc3\", \"Art3\", \"Arvcf\", \"Arvcf\", \"Arxes2\", \"Tns3\", \"Asb8\", \"Aspa\", \"Asrgl1\", \"Atad2\", \"Asxl3\", \"Tnik\", \"Atad2\", \"Tnfrsf21\", \"Atg3\", \"Tnfaip6\", \"Tmpo\", \"Atp1b1\", \"Tmod2\", \"Tmem9b\", \"Tmem63b\", \"Atp8a1\", \"Aurkb\", \"Tmem245\", \"B4galt6\", \"Bach1\", \"Bambi\", \"Tmem117\", \"Tm7sf3\", \"Bcas1\", \"Tk1\", \"Tjp1\", \"Bcl7a\", \"Thtpa\", \"Bicd1\", \"Bin3\", \"Birc5\", \"Tgfbr3\", \"Borcs5\", \"Bpnt1\", \"Tex30\", \"Tdrkh\", \"Tead1\", \"Bzw2\", \"Tcp11l1\", \"Tcf19\", \"Tbx2\", \"Cab39l\", \"Cadm1\", \"Cadm4\", \"Camk2d\", \"Tbc1d10a\", \"Tax1bp3\", \"Capg\", \"Tanc2\", \"Capn5\", \"Taldo1\", \"Tagln2\", \"Taf13\", \"Tacc1\", \"Syt11\", \"Syngr1\", \"Cbx6\", \"Svip\", \"Ccdc13\", \"Ccdc28b\", \"Stx7\", \"Strbp\", \"Ccna2\", \"Ccnd1\", \"Ccnb2\", \"Stmn2\", \"Ccser2\"], \"folder_id\": null, \"folder_parent_id\": null, \"folder_label\": null, \"user_name\": \"Test Armstrong\", \"gene_count\": 159, \"organism\": \"Mus musculus\", \"is_owner\": false}",
                "{\"id\": 202, \"user_id\": 2, \"gctype\": \"unweighted-list\", \"label\": \"demo gene collection\", \"organism_id\": 1, \"ldesc\": null, \"share_id\": \"2fe67c37\", \"is_public\": 1, \"is_domain\": null, \"date_added\": \"2022-02-16 14:27:24\", \"genes\": [\"Pde1c\", \"Nav1\", \"Zfp37\", \"Chl1\", \"Mak16\", \"Gab1\", \"Stard9\", \"Farp2\", \"Galnt6\", \"Hjurp\", \"Sept9\", \"Gm14150\", \"Gm16139\", \"Gm35013\"], \"folder_id\": null, \"folder_parent_id\": null, \"folder_label\": null, \"user_name\": \"Test Armstrong\", \"gene_count\": 14, \"organism\": \"Mus musculus\", \"is_owner\": false}",
            ],
            "pagination": {
                "total_results": 2,
                "current_page": "1",
                "limit": "20",
                "total_pages": 1,
                "next_page": null,
                "prev_page": null
            }
        }
        await route.fulfill({ json });
    });
}

/**
 * Mocks the behavior of searching for gene collections when no results are found.
 * @param {Page} page - The page object used for mocking the route.
 * @returns {Promise<void>}
 */
const mockSearchGeneCollectionsNothingFound = async (page) => {
    await page.route(`${gearBase}/cgi/search_gene_carts.cgi`, async route => {
        const json = {
            "success": 1,
            "problem": "",
            "gene_carts": [],
            "pagination": {
                "total_results": 0,
                "current_page": "1",
                "limit": "20",
                "total_pages": 0,
                "next_page": null,
                "prev_page": null
            }
        }
        await route.fulfill({ json });
    });
}

// ? should i just make this a handler to edit the JSON from mockSearchGeneCollectionsNothingFound before fulfilling?
/**
 * Mocks the behavior of saving a new gene collection.
 * @param {Page} page - The page object for testing.
 * @returns {Promise<void>}
 */
const mockSaveNewGeneCollection = async (page) => {
    await page.route(`${gearBase}/cgi/save_new_genecart_form.cgi`, async route => {
        const json = {
            "id": 2
        };
        await route.fulfill({ json });
    });

    await page.route(`${gearBase}/cgi/search_gene_carts.cgi`, async route => {
        const json = {
            "success": 1,
            "problem": "",
            "gene_carts": [
                JSON.stringify({"id": 2, "user_id": 1, "gctype": "unweighted-list", "label": "new test", "organism_id": 1, "ldesc": null, "share_id": "abcdefg", "is_public": 0, "is_domain": null, "date_added": "2022-02-16 14:27:24", "genes": ["Pou4f3", "Sox2"], "folder_id": null, "folder_parent_id": null, "folder_label": null, "user_name": "Test Armstrong", "gene_count": 14, "organism": "Mus musculus", "is_owner": true}),
            ],
            "pagination": {
                "total_results": 1,
                "current_page": "1",
                "limit": "20",
                "total_pages": 1,
                "next_page": null,
                "prev_page": null
            }
        }
        await route.fulfill({ json });
    })
}

const mockSaveGeneCollectionChanges = async (page) => {
    await page.route(`${gearBase}/cgi/save_genecart_changes.cgi`, async route => {
        const json = {
            "gene_cart": "{\"id\": 334, \"user_id\": 1, \"gctype\": \"unweighted-list\", \"label\": \"updated test\", \"organism_id\": 1, \"ldesc\": \"This is the unweighted description\", \"share_id\": \"06c1cdd3\", \"is_public\": 1, \"is_domain\": null, \"date_added\": \"2023-11-01 17:10:36\", \"genes\": [\"0610040J01Rik\", \"Zkscan1\", \"1500011K16Rik\", \"1700025G04Rik\", \"Zfp60\", \"1810037I17Rik\", \"2010111I01Rik\", \"2300009A05Rik\", \"Zfp24\", \"Zdhhc2\", \"4931406P16Rik\", \"Ypel2\", \"5031439G07Rik\", \"Wwp1\", \"Wsb2\", \"Wls\", \"Wasf1\", \"Vwc2\", \"Vwa1\", \"Aatk\", \"Vps36\", \"Vps37b\", \"Abca8a\", \"Vldlr\", \"Abhd4\", \"Vcl\", \"Utrn\", \"Acsbg1\", \"Usp1\", \"Adam10\", \"Adam17\", \"Adam23\", \"Adamts17\", \"Adamts5\", \"Ung\", \"Unc119\", \"Adgrg6\", \"Ucp2\", \"Ado\", \"Ado\", \"Uchl1\", \"Adra2c\", \"Afap1l2\", \"Ubl3\", \"Ube2e2\", \"Ube2c\", \"Ahr\", \"Ak3\", \"Ak6\", \"Tyms\", \"Akr1b8\", \"Twf1\", \"Alcam\", \"Alad\", \"Tubg1\", \"Ttyh1\", \"Ttyh2\", \"Tst\", \"Tspo\", \"Angpt2\", \"Ank3\", \"Tspan18\", \"Tspan13\", \"Tspan15\", \"Tspan17\", \"Tsc22d4\", \"Anxa2\", \"Trim2\", \"Trim13\", \"Arfip1\", \"Trappc2\", \"Arhgap19\", \"Arhgap24\", \"Arhgap32\", \"Arhgap39\", \"Arhgef10\", \"Arhgef28\", \"Tpst1\", \"Tppp3\", \"Tpp1\", \"Arid5b\", \"Tpd52\", \"Tox\", \"Arl8a\", \"Arnt2\", \"Top2a\", \"Arrdc3\", \"Art3\", \"Arvcf\", \"Arvcf\", \"Arxes2\", \"Tns3\", \"Asb8\", \"Aspa\", \"Asrgl1\", \"Atad2\", \"Asxl3\", \"Tnik\", \"Atad2\", \"Tnfrsf21\", \"Atg3\", \"Tnfaip6\", \"Tmpo\", \"Atp1b1\", \"Tmod2\", \"Tmem9b\", \"Tmem63b\", \"Atp8a1\", \"Aurkb\", \"Tmem245\", \"B4galt6\", \"Bach1\", \"Bambi\", \"Tmem117\", \"Tm7sf3\", \"Bcas1\", \"Tk1\", \"Tjp1\", \"Bcl7a\", \"Thtpa\", \"Bicd1\", \"Bin3\", \"Birc5\", \"Tgfbr3\", \"Borcs5\", \"Bpnt1\", \"Tex30\", \"Tdrkh\", \"Tead1\", \"Bzw2\", \"Tcp11l1\", \"Tcf19\", \"Tbx2\", \"Cab39l\", \"Cadm1\", \"Cadm4\", \"Camk2d\", \"Tbc1d10a\", \"Tax1bp3\", \"Capg\", \"Tanc2\", \"Capn5\", \"Taldo1\", \"Tagln2\", \"Taf13\", \"Tacc1\", \"Syt11\", \"Syngr1\", \"Cbx6\", \"Svip\", \"Ccdc13\", \"Ccdc28b\", \"Stx7\", \"Strbp\", \"Ccna2\", \"Ccnd1\", \"Ccnb2\", \"Stmn2\", \"Ccser2\"], \"folder_id\": null, \"folder_parent_id\": null, \"folder_label\": null, \"user_name\": \"Test Armstrong\", \"gene_count\": 159, \"organism\": \"Mus musculus\", \"is_owner\": false}",
            "success": 1
        };
        await route.fulfill({ json });
    });
}

/**
 * Mocks the delete gene collection functionality.
 * @param {Page} page - The page object.
 * @returns {Promise<void>}
 */
const mockDeleteGeneCollection = async (page) => {
    await page.route(`${gearBase}/cgi/remove_gene_cart.cgi`, async route => {
        const json = {
            "success": 1
        };
        await route.fulfill({ json });
    });

    await page.route(`${gearBase}/cgi/search_gene_carts.cgi`, async route => {
        const json = {
            "success": 1,
            "problem": "",
            "gene_carts": [
                "{\"id\": 332, \"user_id\": 1, \"gctype\": \"weighted-list\", \"label\": \"Weighted test\", \"organism_id\": 1, \"ldesc\": \"This is the weighted description\", \"share_id\": \"b5102349\", \"is_public\": 0, \"is_domain\": null, \"date_added\": \"2023-10-24 14:32:04\", \"genes\": [\"0610025J13Rik\", \"0610030E20Rik\", \"0610030E20Rik\", \"0610033M10Rik\", \"0610040B10Rik\", \"1110002L01Rik\", \"1110008L16Rik\", \"1110008P14Rik\", \"1110017D15Rik\", \"1110020A21Rik\", \"1110032F04Rik\"], \"folder_id\": null, \"folder_parent_id\": null, \"folder_label\": null, \"user_name\": \"Test Armstrong\", \"gene_count\": 11, \"organism\": \"Mus musculus\", \"is_owner\": true}",
                "{\"id\": 202, \"user_id\": 2, \"gctype\": \"unweighted-list\", \"label\": \"demo gene collection\", \"organism_id\": 1, \"ldesc\": null, \"share_id\": \"2fe67c37\", \"is_public\": 1, \"is_domain\": null, \"date_added\": \"2022-02-16 14:27:24\", \"genes\": [\"Pde1c\", \"Nav1\", \"Zfp37\", \"Chl1\", \"Mak16\", \"Gab1\", \"Stard9\", \"Farp2\", \"Galnt6\", \"Hjurp\", \"Sept9\", \"Gm14150\", \"Gm16139\", \"Gm35013\"], \"folder_id\": null, \"folder_parent_id\": null, \"folder_label\": null, \"user_name\": \"Test Armstrong\", \"gene_count\": 14, \"organism\": \"Mus musculus\", \"is_owner\": false}",
            ],
            "pagination": {
                "total_results": 2,
                "current_page": "1",
                "limit": "20",
                "total_pages": 1,
                "next_page": null,
                "prev_page": null
            }
        }
        await route.fulfill({ json });
    })
}

/**
 * Mocks the unweighted gene collection preview.
 *
 * @param {Page} page - The page object.
 * @returns {Promise<void>} - A promise that resolves when the gene collection preview is mocked.
 */
const mockUnweightedGeneCollectionPreview = async (page) => {
    await page.route(`${gearBase}/cgi/get_unweighted_gene_cart_preview.cgi`, async route => {
        const json = {
                "success": 1,
                "gene_info": {
                    "ENSMUSG00000060512": {
                        "gene_symbol": "0610040J01Rik",
                        "product": "RIKEN cDNA 0610040J01 gene"
                    },
                    "ENSMUSG00000029729": {
                        "gene_symbol": "Zkscan1",
                        "product": "zinc finger with KRAB and SCAN domains 1"
                    },
                    "ENSMUSG00000051319": {
                        "gene_symbol": "1500011K16Rik",
                        "product": "RIKEN cDNA 1500011K16 gene"
                    },
                    "ENSMUSG00000032666": {
                        "gene_symbol": "1700025G04Rik",
                        "product": "RIKEN cDNA 1700025G04 gene"
                    },
                    "ENSMUSG00000037640": {
                        "gene_symbol": "Zfp60",
                        "product": "zinc finger protein 60"
                    },
                }
            }
        await route.fulfill({ json });
    })
}

/**
 * Mocks the weighted gene collection preview by intercepting a network request and returning a predefined JSON response.
 * @param {Page} page - The page object representing the browser page.
 * @returns {Promise<void>} - A promise that resolves when the network request is intercepted and fulfilled.
 */
const mockWeightedGeneCollectionPreview = async (page) => {
    await page.route(`${gearBase}/cgi/get_weighted_gene_cart_preview.cgi`, async route => {
        const json = {"preview_json": [], "success": 1, "num_genes": 5486, "weights": ["FC"]}
        await route.fulfill({ json });
    })
}

const mockDownloadUnweightedGeneCollectionMembers = async (page) => {
    await page.route(`${gearBase}/cgi/get_gene_cart_members.cgi`, async route => {
        const json = {
            "gene_symbols": [
                {
                    "id": 111273,
                    "label": "Acbd7"
                },
                {
                    "id": 111274,
                    "label": "Calb1"
                },
                {
                    "id": 111275,
                    "label": "Calm1"
                },
                {
                    "id": 111276,
                    "label": "Calm2"
                },
                {
                    "id": 111277,
                    "label": "Cib2"
                },
                {
                    "id": 111278,
                    "label": "Espn"
                },
                {
                    "id": 111279,
                    "label": "Evl"
                },
                {
                    "id": 111280,
                    "label": "Myo6"
                },
                {
                    "id": 111281,
                    "label": "Pcp4"
                },
                {
                    "id": 111282,
                    "label": "Pou4f3"
                },
                {
                    "id": 111283,
                    "label": "Rasd2"
                },
                {
                    "id": 111284,
                    "label": "Smpx"
                },
                {
                    "id": 111285,
                    "label": "Stard10"
                },
                {
                    "id": 111286,
                    "label": "Tpm1"
                }
            ],
            "success": 1
        }
        await route.fulfill({ json });
    })
}

describe('Gene Collection Manager', function () {
    this.retries(3);
    this.timeout(10000);    // default is 2000

    let browserIndex = 0;

    beforeEach("Setting up", async () => {

        page = await setupBrowserContext(browserIndex);
        await page.goto(gearUrl);
    });

    afterEach("Closing", async () => {
        await teardownBrowserContext();
    });

    for (const b of browsers) {

        describe(`Browser: ${b.name}`, () => {

            after("Switching browsers", () => {
                // Runs once after all tests in this block
                browserIndex++;
            });

            it('Checks that page is returned', () => {
                // Ensure response has table element
                const controlFacet = page.locator("css=#controls_ownership");
                expect(controlFacet).toBeDefined();
            });

            it("should return incorrect password if login failure", async () => {
                // Mock login
                await loginFailure(page, gearUrl);

                // Check that the error message is displayed
                const incorrectPw = page.getByText("Incorrect password");
                await expect(incorrectPw).toBeVisible();
            });

            describe("New Gene Collection", () => {
                it("should not appear if not logged in", async () => {
                    // Check that the form is hidden
                    await expect(page.getByRole("button", {name: "Create new gene collection"})).not.toBeVisible();
                });

                describe("logged in", () => {
                    beforeEach("logging in", async () => {
                        await login(page, gearUrl);
                    });

                    it('should validate #new_collection_label input on blur', async function () {
                        await page.getByRole("button", {name: "Create new gene collection"}).click();
                        await page.getByRole("button", {name: "Paste list of genes"}).click();

                        // Trigger the blur event
                        const label = await page.locator("css=#new_collection_label");
                        await label.click();
                        await page.locator("css=#new_collection_ldesc").click();

                        // Check that the correct class was added/removed from #new_collection_label
                        // Check that the helper text was added/removed
                        await expect(label).toHaveClass(/is-danger/);
                        await expect(page.getByText("Please enter a value")).toBeVisible();
                    });

                    it('should validate #new_collection_label input on save', async function () {
                        await page.getByRole("button", {name: "Create new gene collection"}).click();
                        await page.getByRole("button", {name: "Paste list of genes"}).click();

                        // Trigger the click event
                        await page.locator("css=#btn_new_collection_save").click();

                        await expect(page.locator("css=#new_collection_label")).toHaveClass(/is-danger/);

                        const parent = page.locator("css=.control:has(#new_collection_label)")
                            .filter({ has: page.getByText("Please enter a value") });

                        await expect(parent).toBeVisible();
                    });

                    it("should validate #new_collection_pasted_genes input on save", async function () {
                        await page.getByRole("button", {name: "Create new gene collection"}).click();
                        await page.getByRole("button", {name: "Paste list of genes"}).click();

                        // Trigger the click event
                        await page.locator("css=#btn_new_collection_save").click();

                        await expect(page.locator("css=#new_collection_label")).toHaveClass(/is-danger/);

                        const parent = page.locator("css=.control:has(#new_collection_pasted_genes)")
                            .filter({ has: page.getByText("Please enter a value") });

                        await expect(parent).toBeVisible();
                    });

                    it("should validate #new_collection_file input on save", async function () {
                        await page.getByRole("button", {name: "Create new gene collection"}).click();
                        await page.locator("css=#btn_gc_upload_unweighted_list").click();

                        // Trigger the click event
                        await page.locator("css=#btn_new_collection_save").click();

                        await expect(page.locator("css=.file:has(#new_collection_file)")).toHaveClass(/is-danger/);
                        await expect(page.getByText("Please select a file")).toBeVisible();
                    });

                    it('should validate organism select on save', async function () {
                        await page.getByRole("button", {name: "Create new gene collection"}).click();
                        await page.getByRole("button", {name: "Paste list of genes"}).click();

                        // Trigger the click event
                        await page.locator("css=#btn_new_collection_save").click();

                        await expect(page.locator("css=.select:has(#new_collection_organism_id)")).toHaveClass(/is-danger/);
                        await expect(page.getByText("Please select an organism")).toBeVisible();
                    });

                    it("should hide #new_collection_form_c on cancel", async function () {
                        await page.getByRole("button", {name: "Create new gene collection"}).click();
                        await page.getByRole("button", {name: "Paste list of genes"}).click();

                        // Trigger the click event
                        await page.locator("css=#btn_new_collection_cancel").click();

                        // Check that the form is hidden
                        await expect(page.locator("css=#new_collection_form_c")).not.toBeVisible();
                    });

                    it("should show new gene collection in results after save", async function () {
                        await mockSaveNewGeneCollection(page);

                        await page.getByRole("button", {name: "Create new gene collection"}).click();
                        await page.getByRole("button", {name: "Paste list of genes"}).click();

                        // Populate label, genes, and organism
                        await page.locator("css=#new_collection_label").fill("new test");
                        await page.locator("css=#new_collection_pasted_genes").fill("Pou4f3, Sox2");
                        await page.locator("css=#new_collection_organism_id").selectOption({ label:"Mouse" });

                        // Trigger the click event
                        await page.locator("css=#btn_new_collection_save").click();

                        // Check that the form is hidden
                        await expect(page.locator("css=#new_collection_form_c")).not.toBeVisible();

                        // Check that the search results were updated
                        await expect(page.locator("css=#gc_count_label")).toContainText("result");  // 1 result

                        // Check that the new gene collection is in the results
                        const resultsTable = page.locator("css=#results_table");
                        await expect(resultsTable.getByText("new test")).toBeVisible();
                    });

                });

            });

            describe("Gene Collection Search", () => {
                describe("Search Results", () => {

                    beforeEach("Mocking search", async () => {
                        await mockSearchGeneCollections(page);
                    });

                    it('should clear search terms on #search_clear click', async () => {
                        // Set the #search_terms value
                        await page.locator("css=#search_terms").fill("test");
                        await page.keyboard.up("t")
                        // Trigger the click event
                        await page.locator("css=#search_clear").click();
                        // Check that the #search_terms value was cleared
                        await expect(page.locator("css=#search_terms")).toHaveValue("");
                    });

                    it('should search after enter clicked for #search_terms', async () => {
                        // Trigger the keyup event
                        await page.locator("css=#search_terms").type("test");
                        await page.keyboard.up("Enter");
                        // Check that submitSearch() was called if the key was "Enter"
                        await expect(page.locator("css=#gc_count_label")).toContainText("results");
                    });

                    it('should update search results on #sort_by change', async () => {
                        // Perform a search
                        await page.locator("css=#search_terms").fill("test");
                        await page.keyboard.up("Enter");
                        // Trigger the change event
                        await page.locator("css=#sort_by").selectOption("Gene cart title");
                        // Check that the search results were updated
                        await expect(page.locator("css=#gc_count_label")).toContainText("results");
                    });

                    it("should deselect other #controls_date_added facets if one is selected", async () => {
                        const anyTimeOption = page.getByText("Any time");
                        const lastWeekOption = page.getByText("Within last week");
                        const lastMonthOption = page.getByText("Within last month");

                        // Selecting should deselect "All"
                        await lastWeekOption.click();
                        await expect(lastWeekOption).toHaveClass(/js-selected/);
                        await expect(anyTimeOption).not.toHaveClass(/js-selected/);

                        // Switching should deselect the previous one
                        await lastMonthOption.click();
                        await expect(lastMonthOption).toHaveClass(/js-selected/);
                        await expect(lastWeekOption).not.toHaveClass(/js-selected/);
                    });

                    it("should deselect non-All facets if All is selected", async () => {
                        await login(page, gearUrl);

                        const groupAffiliatedOption = page.getByText("Group-affiliated collections");
                        const yourCollectionsOption = page.getByText("Your collections");
                        const allFacet = page.locator("css=#controls_ownership .js-all-selector");

                        // Selecting should deselect "All"
                        await groupAffiliatedOption.click();
                        // TODO: detecting the class seems to be flaky when all tests are run
                        await expect(groupAffiliatedOption).toHaveClass(/js-selected/);
                        await expect(yourCollectionsOption).toHaveClass(/js-selected/);
                        await expect(allFacet).not.toHaveClass(/js-selected/);

                        // Switching to All should deselect all others
                        await allFacet.click();
                        await expect(allFacet).toHaveClass(/js-selected/);
                        await expect(groupAffiliatedOption).not.toHaveClass(/js-selected/);
                        await expect(yourCollectionsOption).not.toHaveClass(/js-selected/);
                    });

                    it("should only show public gene collections when not logged in", async () => {
                        await mockSearchGeneCollectionsNotLoggedIn(page);
                        // Only 2 public gene collections from our mocked response
                        await expect(page.locator("css=.js-gc-list-element")).toHaveCount(2);
                    })

                    describe("Gene collection actions", () => {

                        it("should show the table view when that button is clicked", async () => {
                            await page.locator("css=#btn_table_view").click();
                            await expect(page.locator("css=#results_table")).toBeVisible();
                            await expect(page.locator("css=#results_list_div")).not.toBeVisible();
                            await expect(page.getByText("This is the unweighted description")).not.toBeVisible();
                            await expect(page.getByText("This is the weighted description")).not.toBeVisible();
                        });

                        it("should expand/collapse all results when those buttons are clicked", async () => {
                            // Expand all
                            await page.locator("css=#btn_list_view_expanded").click();
                            await expect(page.getByText("This is the unweighted description")).toBeVisible();
                            await expect(page.getByText("This is the weighted description")).toBeVisible();
                            await expect(page.locator("css=#results_table")).not.toBeVisible();
                            // Collapse all
                            await page.locator("css=#btn_list_view_compact").click();
                            await expect(page.getByText("This is the unweighted description")).not.toBeVisible();
                            await expect(page.getByText("This is the weighted description")).not.toBeVisible();
                            await expect(page.locator("css=#results_table")).not.toBeVisible();

                        });

                        it("should expand/collapse individual results when those buttons are clicked", async () => {
                            // switch to list view
                            await page.locator("css=#btn_list_view_compact").click();

                            // Expand
                            await page.locator("css=#result_gc_id_334 .js-expand-box").click();
                            await expect(page.getByText("This is the unweighted description")).toBeVisible();
                            await expect(page.getByText("This is the weighted description")).not.toBeVisible();

                            // Collapse
                            await page.locator("css=#result_gc_id_334 .js-expand-box").click();
                            await expect(page.getByText("This is the unweighted description")).not.toBeVisible();
                        })

                        it("edit and delete buttons should not appear if not logged in", async () => {
                            await mockSearchGeneCollectionsNotLoggedIn(page);
                            // Check that the form is hidden
                            const resultsTable = page.locator("css=#results_table");
                            await expect(resultsTable.getByText("Unweighted test")).toBeVisible();
                            await expect(page.locator("css=#result_gc_id_334 .js-edit-gc")).not.toBeVisible();
                            await expect(page.locator("css=#result_gc_id_334 .js-delete-gc")).not.toBeVisible();
                        });

                        it("should show unweighted genes in table when preview genes button is clicked", async () => {
                            // switch to list view
                            await page.locator("css=#btn_list_view_compact").click();

                            await mockUnweightedGeneCollectionPreview(page);

                            const previewBtn = page.locator("css=#result_gc_id_334 .js-preview-genes-button-container");
                            const previewContainer = page.locator("css=#result_gc_id_334 .js-preview-genes-container");

                            await expect(previewBtn).toHaveText("159 genes");
                            await previewBtn.getByText("159 genes").click();
                            await expect(previewBtn).toHaveText("Hide");

                            await expect(previewContainer.getByText("0610040J01Rik")).toBeVisible();
                            await expect(previewContainer.getByText("ENSMUSG00000051319")).toBeVisible();   // second gene in table
                            await expect(previewContainer.getByText("RIKEN cDNA 1700025G04 gene")).toBeVisible();   // third gene in table

                            await previewBtn.getByText("Hide").click();
                            await expect(previewBtn).toHaveText("159 genes");
                        })

                        it("should show weighted gene infomation when preview genes button is clicked", async () => {
                            // switch to list view
                            await page.locator("css=#btn_list_view_compact").click();

                            await mockWeightedGeneCollectionPreview(page);

                            const previewBtn = page.locator("css=#result_gc_id_332 .js-preview-genes-button-container");
                            const previewContainer = page.locator("css=#result_gc_id_332 .js-preview-genes-container");

                            await expect(previewBtn).toHaveText("Info");
                            await previewBtn.getByText("Info").click();
                            await expect(previewBtn).toHaveText("Hide");
                            await expect(previewContainer.getByText("5486")).toBeVisible(); // genes
                            await expect(previewContainer.getByText("1")).toBeVisible();   // num weights
                            await expect(previewContainer.getByText("FC")).toBeVisible();   // weight labels

                            await previewBtn.getByText("Hide").click();
                            await expect(previewBtn).toHaveText("Info");
                        })

                        describe("login required", () => {
                            beforeEach("logging in", async () => {
                                await login(page, gearUrl);

                                // switch to list view
                                await page.locator("css=#btn_list_view_compact").click();
                            });

                            it("edit and delete buttons should not appear for datasets user does not own", async () => {
                                // Check that the form is hidden
                                const resultsList = page.locator("css=#results_list_div");
                                await expect(resultsList.getByText("Unweighted test")).toBeVisible();
                                await expect(page.locator("css=#result_gc_id_202 .js-edit-gc")).not.toBeVisible();
                                await expect(page.locator("css=#result_gc_id_202 .js-delete-gc")).not.toBeVisible();
                            });

                            it("should show edit form when edit button is clicked", async () => {
                                await page.locator("css=#result_gc_id_334 .js-edit-gc").click();
                                await expect(page.locator("css=#result_gc_id_334_editable_title")).toBeVisible();
                            });

                            it("should update gene collection when save button is clicked", async () => {
                                await page.locator("css=#result_gc_id_334 .js-edit-gc").click();
                                await page.locator("css=#result_gc_id_334_editable_title").fill("updated test");

                                await mockSaveGeneCollectionChanges(page);

                                await page.locator("css=#result_gc_id_334 .js-edit-gc-save").click();
                                await expect(page.locator("css=#result_gc_id_334_editable_title")).not.toBeVisible();
                                // Check that the search results were updated
                                await expect(page.locator("css=#gc_count_label")).toContainText("results");
                                // Check that the updated gene collection is in the results
                                const resultsList = page.locator("css=#results_list_div");
                                await expect(resultsList.getByText("updated test")).toBeVisible();
                            });

                            it("should hide edit form when cancel button is clicked", async () => {
                                await page.locator("css=#result_gc_id_334 .js-edit-gc").click();
                                await page.locator("css=#result_gc_id_334 .js-edit-gc-cancel").click();
                                await expect(page.locator("css=#result_gc_id_334_editable_title")).not.toBeVisible();
                            });

                            it("should delete gene collection when delete button is clicked", async () => {
                                // deleting gene collection with ID 334
                                await page.locator("css=#result_gc_id_334 .js-delete-gc").click();
                                await expect(page.getByText("Remove collection")).toBeVisible();

                                await mockDeleteGeneCollection(page);

                                await page.locator("css=#confirm_gc_delete").click();
                                await expect(page.getByText("Remove collection")).not.toBeVisible();
                                // Check that the search results were updated and the deleted gene collection is gone
                                await expect(page.locator("css=#gc_count_label")).toContainText("results");
                                await expect(page.locator("css=.js-gc-list-element")).toHaveCount(2);
                                await expect(page.getByText("Unweighted test")).not.toBeVisible();
                            });

                            it("should cancel delete gene collection when cancel button is clicked", async () => {
                                await page.locator("css=#result_gc_id_334 .js-delete-gc").click();
                                await expect(page.getByText("Remove collection")).toBeVisible();
                                await page.locator("css=#cancel_gc_delete").click();
                                await expect(page.getByText("Remove collection")).not.toBeVisible();
                                // Ensure the gene collection is still there
                                const resultsList = page.locator("css=#results_list_div");
                                await expect(resultsList.getByText("Unweighted test")).toBeVisible();
                            });

                            // TODO: How to test copied to clipboard?
                            it.skip("should copy share link when copy button is clicked", async () => {
                                await page.locator("css=#result_gc_id_334 .js-share-gc").click();
                                await expect(page.getByText("URL copied to clipboard")).toBeVisible();
                            });

                            it("should redirect to gene search view when view button is clicked", async () => {
                                const viewBtn = page.locator("css=#result_gc_id_334 .js-view-gc");
                                const shareId = await viewBtn.getAttribute("value");
                                await viewBtn.click();
                                // On new page now.  Share ID should be in the URL
                                expect(page.url()).toContain(shareId);
                            });

                            it("should download unweighted gene collection when download button is clicked", async () => {
                                await mockDownloadUnweightedGeneCollectionMembers(page);
                                // Start waiting for download before clicking. Note no await.
                                const downloadPromise = page.waitForEvent('download');

                                await page.locator("css=#result_gc_id_334 .js-download-gc").click();
                                const download = await downloadPromise;
                                // Wait for the download process to complete and save the downloaded file somewhere.
                                await download.saveAs(`/tmp/${download.suggestedFilename()}`);
                                // TODO: How to test downloaded file contents?
                            })

                            it("should download weighted gene collection when download button is clicked", async () => {
                                // don't need to mock here has it outputs an attachment

                                // Start waiting for download before clicking. Note no await.
                                const downloadPromise = page.waitForEvent('download');

                                await page.locator("css=#result_gc_id_332 .js-download-gc").click();
                                const download = await downloadPromise;
                                // Wait for the download process to complete and save the downloaded file somewhere.
                                await download.saveAs(`/tmp/${download.suggestedFilename()}`);
                            })
                        });

                    });

                });

                describe("No Search Results", () => {
                    beforeEach("Mocking search", async () => {
                        await mockSearchGeneCollectionsNothingFound(page);
                    });

                    it('should display "No results found" when no results are found', async () => {
                        // Perform a search
                        await page.locator("css=#search_terms").fill("test");
                        await page.keyboard.up("Enter");
                        // Check that the search results were updated
                        await expect(page.getByText("No results found")).toBeVisible();
                    });
                });

            });

        });
    };
});