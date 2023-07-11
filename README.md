## cancerHotspots

R data package with mutational hotspots in cancer. The data has been harvested from cancerhotspots.org. Location of each hotspot is provided with HGVSp/HGVSc coordinates. Imporantly, a number of splice site hotspots have been curated with information regarding their relative position with respect to the exon/intron boundary (i.e. HGVSc).

### Installation

`remotes::install_github('sigven/cancerHotspots')`

### Usage

All mutation hotspots - records in _long_ format (one record per hotspot per cancer type):

- `cancerHotspots::cancer_hotspots[['long']]`

All mutation hotspots - records in _wide_ format:

- `cancerHotspots::cancer_hotspots[['wide']]`

Metadata (citation, license etc:)

- `cancerHotspots::cancer_hotspots[['metadata']]`

### Important note

If you use the dataset provided with **cancerHotspots**, make sure you properly cite the original publication:

- [Chang et al., Cancer Discov, 2018](https://pubmed.ncbi.nlm.nih.gov/29247016/)
