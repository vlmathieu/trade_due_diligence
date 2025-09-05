# trade_due_diligence
A repository to analyse the impact of due diligence policies on trade.

## Folder structure
The folder is organized as follows:

```bash
├── README.md
├── resources
│   └── inhouse
│       └── results-survey857139.csv
├── results
└── workflow
    ├── Snakefile
    ├── envs
    ├── rules
    ├── sandbox
    └── scripts
```


README.md provides information on the repository structuration and explains the data analysis workflow.
 
The workflow code goes into a subfolder `workflow`, while the configuration is stored in a subfolder `config`. Inside of the workflow subfolder, the central `Snakefile` marks the entrypoint of the workflow (it will be automatically discovered when running snakemake from the root of above structure). In addition to the central `Snakefile`, rules are stored in a modular way, using the subfolder `workflow/rules`. Such modules should end with `.smk`, the recommended file extension of Snakemake. Further, scripts are stored in a subfolder `workflow/scripts`. Conda environments are stored in the subfolder `workflow/envs` (they are kept as finegrained as possible to improve transparency and maintainability).

All output files generated in the workflow are stored under `results`, unless they are rather retrieved `resources`, in which case they should be stored under resources. The latter subfolder also contains small resources that shall be delivered along with the workflow via git.

# To-do list
- [x] Regarder les HS code à retenir :
    - [x] Ceux couvert par les politiques (baseline EUTR)
    - [x] Ceux avec une correspondance FAO pour estimatio, flux i->i
- [x] Une fois fait, réfléchir à la méthode d'aggrégation des productions FAO
- [] Check if due diligence policies are the only policy tools implemented to fight illegal trade of wood
- [] Check the share of deforestation and forest degradation due to illegal activities (and if possible due to illegal timber trade)

# Ideas for discussion

*Was it worth it?*

Illegal timber trade is an important revenue from environmental crime (half of illicit proceed from environmental crime), but most proably not the main driver of deforestation.

Deforestation results from structural domestic dynamics, that can be legal or illegal: commercial agriculture, pasture, shifting cultivation, other subsistence agriculture, managed forests, wildfire, natural disturbances, roads, mining, commercial oil palm.

Dummett and Blundell, 2021: Almost two-thirds (60 percent) of tropical deforestation between 2013 and 2019 was driven by commercial agriculture. At least 69 percent of agro-conversion was conducted in violation of national laws and regulations, and this is likely an underestimate. More than 31 percent of agricultural commodities linked to deforestation were exported, raising significant concerns about their association with illegal deforestation. https://www.forest-trends.org/wp-content/uploads/2021/05/Illicit-Harvest-Complicit-Goods_rev.pdf

From the illegal part, a part is converted to trade. From Dummett and Blundell (2021): 0.6 x 0.69 x X 0.31 = 0.13 ~ 13%

Focusing on wood products also shows limits: most of deforestation is linked to agriculture.

Due diligence policies: high cost of implementation, high opportunity costs, non-tariff barrier for everyone, so compliance cost for everyone + undifferentiated effect whenever a country is risky/dirty or not + redistribution effect.

At a time when occidental countries need biomass to achieve sustainable goals and are dependant on trade, did it worth it?
At least, it reduces the responsibility of occidental countries.

Speak about workaround strategy (as perspective):
- Mix wood products that are regulated with other product so that they change HS codes and are not regulated anymore
- Deviation of trade flows to more relaxed country
- More domestic consumption
- More processed wood products to blur traceability and pass regulation easily
