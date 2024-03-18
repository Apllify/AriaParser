import pandas as pd


path = "/Users/anjamatic/AriaParser/output_data/hmqc_14/"
methods = ["_aarrs"]
n = len(methods)
total_sdist = 0
num_dist = 0

total_sm = 0
total_sp = 0

num_rhos = 0

results = {"sdist":[], "sm":[], "sp":[], "rhos": []}

for i in range(10, 11):
    for sol in methods:
        df = pd.read_csv(f'{path}{i}{sol}_sdists_hmqnoe_14')
        total_sdist +=df['sdist'].sum()
        num_dist = len(df)

        df = pd.read_csv(f'{path}{i}{sol}_sm_hmqnoe_14')
        total_sm +=df['sm'].sum()
        num_rhos = len(df)

        df = pd.read_csv(f'{path}{i}{sol}_sp_hmqnoe_14')
        total_sp +=df['sp'].sum()

    results["sdist"].append(total_sdist/(n*num_dist))
    results["sp"].append( total_sp/((n*num_rhos)))
    results["sm"].append(total_sm/((n*num_rhos)))
    results["rhos"].append(num_rhos)
pd.DataFrame(results).to_csv("/Users/anjamatic/AriaParser/means")   


