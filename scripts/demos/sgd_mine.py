#Automatically checks for gene-drug associations on SGD


from intermine.webservice import Service
service = Service("http://yeastmine.yeastgenome.org/yeastmine/service")

drugs = ['miconazole',
'tamoxifen',
'benomyl',
'fluconazole',
'ketoconazole',
'camptothecin',
'bisantrene',
'cycloheximide',
'itraconazole',
'beauvericin']

for drug in drugs:
    print drug
    query = service.new_query("Phenotype")
    query.add_view(
        "genes.primaryIdentifier", "genes.secondaryIdentifier", "genes.symbol",
        "genes.qualifier", "genes.sgdAlias", "experimentType", "mutantType",
        "observable", "qualifier", "allele", "strainBackground", "chemical",
        "condition", "details", "reporter", "publications.pubMedId",
        "publications.citation"
    )
    query.add_sort_order("Phenotype.experimentType", "ASC")
    query.add_constraint("observable", "=", "Resistance to chemicals", code = "A")
    query.add_constraint("qualifier", "=", "decreased", code = "B")
    query.add_constraint("chemical", "CONTAINS", drug, code = "C")


    results = []
    for row in query.rows():
        results.append(row["genes.symbol"])#, row["genes.secondaryIdentifier"], row["genes.symbol"], \
            #row["genes.qualifier"], row["genes.sgdAlias"], row["experimentType"], row["mutantType"], \
            #row["observable"], row["qualifier"], row["allele"], row["strainBackground"], \
            #row["chemical"], row["condition"], row["details"], row["reporter"], \
            #row["publications.pubMedId"], row["publications.citation"])

    #print results


    gene_list = ['SNQ2','PDR5','YBT1','YCF1','YOR1']
    for gene in gene_list:
        print gene, gene in results