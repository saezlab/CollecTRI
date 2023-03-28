
if(!require("OmnipathR")) remotes::install_github('saezlab/OmnipathR', upgrade='never')

# load helper functions
source('scripts/helper/get_TFcategories.R')
source('scripts/helper/get_networks.R')

# Load TF categories
get_TFcategory('data/TFcategory')
generate_TFlist('data/TFcategory')

# load networks
get_data('data/networks')
get_data('data/networks', filterTFs = F)
