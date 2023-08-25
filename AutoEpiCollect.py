import pandas as pd
import os
import openpyxl
import subprocess
from io import StringIO
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import time
from time import sleep
from selenium import webdriver
from selenium.webdriver.common.by import By
from selenium.webdriver.support.wait import WebDriverWait
from selenium.webdriver.support import expected_conditions as ec
from selenium.common.exceptions import NoSuchElementException
from selenium.webdriver import ChromeOptions
from sklearn import linear_model


def get_gene_sequence(target_gene):
    options = webdriver.ChromeOptions()
    options.add_argument("--headless=new")
    driver = webdriver.Chrome(f"{parent_dir}/chromedriver_mac_arm64/chromedriver", options=options)
    driver.get('https://www.uniprot.org/')

    sleep(1)

    search_box = driver.find_element_by_xpath('//input[@type="text"]')
    search_box.send_keys(target_gene)

    sleep(0.5)

    submit_button = driver.find_element_by_xpath('//*[@type="submit"]')
    submit_button.click()

    sleep(2)

    try:
        choose_table = driver.find_element_by_xpath('/html/body/form/div/span/label[2]/input')
    except NoSuchElementException:
        print("No format selector")
    else:
        choose_table.click()

    sleep(0.5)

    view_results_button = driver.find_element_by_xpath('/html/body/form/div/section/button')
    view_results_button.click()

    sleep(1)

    filter_human = driver.find_element_by_xpath('/html/body/div[1]/div/div/div/aside/div/ul/li[2]/div/ul/li[1]/a')
    filter_human.click()

    sleep(1)

    gene_button = driver.find_element_by_xpath('/html/body/div[1]/div/div/div/main/table/tbody/tr[1]/td[2]/span/a')
    gene_button.click()

    sleep(2)

    download_button = driver.find_element_by_xpath('/html/body/div[1]/div/div/div/main/div/div[2]/div/div/button')
    download_button.click()

    sleep(0.5)

    fasta_button = driver.find_element_by_xpath('/html/body/div[1]/div/div/div/main/div/div[2]/div/div/div/ul/li[2]/a')
    fasta_button.click()

    sleep(2)

    child_tab = driver.window_handles[1]
    driver.switch_to.window(child_tab)

    result_fasta = driver.find_element_by_xpath('/html/body/pre').text
    with open(f"{target_gene}.fasta", "w") as g:
        g.write(result_fasta)

    for tab in driver.window_handles:
        driver.switch_to.window(tab)
        driver.close()


def make_mutant_genes(mutant_list, gene_seq):
    for m in mutant_list:
        file_out = f"{parent_dir}/mutant_gene_fastas/{m}.fasta"
        os.makedirs(os.path.dirname(file_out), exist_ok=True)
        with open(file_out, "w") as f:
            for seq_record in SeqIO.parse(open(gene_seq, mode='r'), 'fasta'):
                wild_type = seq_record.seq
                new = m[-1]
                loc = m[:-1]
                loc = loc[1:]
                loc = int(loc) - 1
                new_seq = wild_type[:loc] + new + wild_type[loc + 1:]
                mutant_seq = SeqRecord(seq=new_seq, id=seq_record.id, description=seq_record.description)
                # write new fasta file
                r = SeqIO.write(mutant_seq, f, 'fasta')
                if r != 1:
                    print('Error while writing sequence:  ' + seq_record.id)


def get_epitopes_ba(mutant_list, mhc, mut_can_dict):
    if mhc == "I":
        file_out = f"mutant_gene_all_epitopes_mhci.xlsx"
        epitopes_dict = {}
        with pd.ExcelWriter(file_out, engine='openpyxl') as writer:
            for m in mutant_list:
                sequence = ""
                epitope_lengths = ""
                fasta_file = f"{parent_dir}/mutant_gene_fastas/{m}.fasta"
                for seq_record in SeqIO.parse(open(fasta_file, mode='r'), 'fasta'):
                    sequence = str(seq_record.seq)
                    with open("MHCI_HLA_input.txt", "r") as h:
                        alleles = h.read() * 2
                        alleles = alleles[:-1]
                    epitope_lengths = "9," * 27 + "10," * 26 + "10"

                mhc_i = subprocess.run(["curl", "--data",
                                        f'method=netmhcpan_ba&sequence_text={sequence}&allele={alleles}&length={epitope_lengths}',
                                        "http://tools-cluster-interface.iedb.org/tools_api/mhci/"], capture_output=True, text=True)
                output = mhc_i.stdout

                table = StringIO(output)
                df = pd.read_table(table, sep=r"\s+")
                df = df.rename(columns={"ic50": "binding affinity (nM)"})
                column_titles = ["allele", "seq_num", "start", "end", "length", "peptide", "core", "icore", "percentile_rank", "binding affinity (nM)"]
                df = df.reindex(columns=column_titles)
                for i in range(df.shape[0]):
                    affinity = float(df["binding affinity (nM)"][i])
                    if affinity <= 50:
                        df.at[i, "binding score"] = "STRONG"
                    elif 50 < affinity <= 500:
                        df.at[i, "binding score"] = "NORMAL"
                    elif 500 < affinity <= 5000:
                        df.at[i, "binding score"] = "WEAK"
                    else:
                        df.at[i, "binding score"] = "N/A"
                print(df)
                epitopes_dict[m] = df
                print(f"{m} done")
                df.to_excel(writer, header=True, index=False, sheet_name=m)
                pd.DataFrame(mut_can_dict[m], columns=["Cancers"]).to_excel(writer, header=True, index=False, sheet_name=m, startcol=df.shape[1] + 1, startrow=0)
    else:
        file_out = f"mutant_gene_all_epitopes_mhcii.xlsx"
        epitopes_dict = {}
        with pd.ExcelWriter(file_out, engine='openpyxl') as writer:
            for m in mutant_list:
                sequence = ""
                epitope_lengths = ""
                fasta_file = f"{parent_dir}/mutant_gene_fastas/{m}.fasta"
                # fasta_file = f"./mutant_gene_fastas/{m}.fasta"
                for seq_record in SeqIO.parse(open(fasta_file, mode='r'), 'fasta'):
                    sequence = str(seq_record.seq)
                    with open("MHCII_HLA_input.txt", "r") as h:
                        alleles = h.read()
                        alleles = alleles[:-1]
                    # epitope_lengths = "12," * 27 + "13," * 27 + "14," * 27 + "15," * 27 + "16," * 27 + "17," * 27 + "18," * 26 + "18"
                    epitope_lengths = "15," * 26 + "15"

                mhc_ii = subprocess.run(["curl", "--data",
                                        f'method=netmhciipan&sequence_text={sequence}&allele={alleles}&length={epitope_lengths}',
                                        "http://tools-cluster-interface.iedb.org/tools_api/mhcii/"],
                                       capture_output=True, text=True)
                output = mhc_ii.stdout

                table = StringIO(output)
                df = pd.read_table(table, sep=r"\s+")
                df = df.rename(columns={"ic50": "binding affinity (nM)"})
                column_titles = ["allele", "seq_num", "start", "end", "length", "core_peptide", "peptide", "rank",
                                 "binding affinity (nM)"]
                df = df.reindex(columns=column_titles)
                for i in range(df.shape[0]):
                    affinity = float(df["binding affinity (nM)"][i])
                    if affinity <= 50:
                        df.at[i, "binding score"] = "STRONG"
                    elif 50 < affinity <= 500:
                        df.at[i, "binding score"] = "NORMAL"
                    elif 500 < affinity <= 5000:
                        df.at[i, "binding score"] = "WEAK"
                    else:
                        df.at[i, "binding score"] = "N/A"
                print(df)
                epitopes_dict[m] = df
                df.to_excel(writer, header=True, index=False, sheet_name=m)
                pd.DataFrame(mut_can_dict[m], columns=["Cancers"]).to_excel(writer, header=True, index=False, sheet_name=m, startcol=df.shape[1] + 1, startrow=0)
                print(f"{m} done")
    return epitopes_dict


def get_mutant_epitopes(mutant_list, mhc, all_epitopes_dict, mut_can_dict):
    if mhc == "I":
        file_out = f"mutant_gene_mut_epitopes_mhci.xlsx"
        epitopes_dict = {}
        with pd.ExcelWriter(file_out, engine='openpyxl') as writer:
            for m in mutant_list:
                fasta_file = f"{parent_dir}/mutant_gene_fastas/{m}.fasta"
                # df = pd.read_excel(f"mutant_gene_all_epitopes_mhci_{c}.xlsx", sheet_name=m)
                df = all_epitopes_dict[m].copy()
                for seq_record in SeqIO.parse(open(fasta_file, mode='r'), 'fasta'):
                    sequence = seq_record.seq
                    loc = m[:-1]
                    loc = loc[1:]
                    loc = int(loc) - 1
                    target_epitope = sequence[loc - 9:loc] + sequence[loc:loc + 10]
                    print(m + ": " + target_epitope)
                    index = 0
                    bad_epitopes_indexes = []
                    for pep in df["peptide"]:
                        if pep not in target_epitope or df["start"][index] == loc + 2 or df["end"][index] == loc:
                            bad_epitopes_indexes.append(index)
                        index += 1
                    df_dropped = df.drop(bad_epitopes_indexes).reset_index(drop=True)
                    df_dropped.to_excel(writer, header=True, index=False, sheet_name=m)
                    epitopes_dict[m] = df_dropped.reset_index(drop=True)
                    pd.DataFrame(mut_can_dict[m], columns=["Cancers"]).to_excel(writer, header=True, index=False, sheet_name=m, startcol=df_dropped.shape[1] + 1, startrow=0)
    else:
        file_out = f"mutant_gene_mut_epitopes_mhcii.xlsx"
        epitopes_dict = {}
        with pd.ExcelWriter(file_out, engine='openpyxl') as writer:
            for m in mutant_list:
                fasta_file = f"{parent_dir}/mutant_gene_fastas/{m}.fasta"
                # df = pd.read_excel(f"mutant_gene_all_epitopes_mhcii_{c}.xlsx", sheet_name=m)
                df = all_epitopes_dict[m].copy()
                for seq_record in SeqIO.parse(open(fasta_file, mode='r'), 'fasta'):
                    sequence = seq_record.seq
                    loc = m[:-1]
                    loc = loc[1:]
                    loc = int(loc) - 1
                    target_epitope = sequence[loc - 14:loc] + sequence[loc:loc + 15]
                    print(m + ": " + target_epitope)
                    bad_epitopes_indexes = []
                    loc += 1
                    for i in range(df.shape[0]):
                        start_index = int(df["start"][i])
                        end_index = int(df["end"][i])
                        pep = df["peptide"][i]
                        if pep not in target_epitope or not start_index <= loc <= end_index:
                            bad_epitopes_indexes.append(i)
                    df_dropped = df.drop(bad_epitopes_indexes).reset_index(drop=True)
                    df_dropped.to_excel(writer, header=True, index=False, sheet_name=m)
                    epitopes_dict[m] = df_dropped
                    pd.DataFrame(mut_can_dict[m], columns=["Cancers"]).to_excel(writer, header=True, index=False, sheet_name=m, startcol=df_dropped.shape[1] + 1, startrow=0)
    return epitopes_dict


def get_peptides(point_mutants, mut_epitopes_dict, mhc):
    for m in point_mutants:
        file_out = f"{parent_dir}/Sequences/{m}peptides_{mhc}.txt"
        os.makedirs(os.path.dirname(file_out), exist_ok=True)
        df = mut_epitopes_dict[m]
        peptide_set = set()
        for pep in df.loc[:, "peptide"]:
            peptide_set.add(pep)
        with open(file_out, "w") as f_out:
            for pep in peptide_set:
                f_out.write(pep + "\n")
        file_out = f"{parent_dir}/Sequences/{m}peptides_{mhc}.fasta"
        os.makedirs(os.path.dirname(file_out), exist_ok=True)
        with open(file_out, "w") as f_out:
            count = 0
            for pep in peptide_set:
                heading = f">Epitope{count}"
                f_out.write(heading + "\n" + pep + "\n")
                count += 1


def get_local_immunogenicity_mhci():
    completed_run = subprocess.run(["python", immunogenicity_file, peptide_file], capture_output=True, text=True)
    output = completed_run.stdout
    output = output[output.find("peptide"):]
    table = StringIO(output)
    df = pd.read_table(table, sep=",")

    for i in range(df.shape[0]):
        result = df["score"][i]
        e = df["peptide"][i]
        current_df.loc[current_df["peptide"] == e, "immunogenicity"] = float(result)


def get_immunogenicity_mhcii(peptide_list):
    options = webdriver.ChromeOptions()
    options.add_argument("--headless=new")
    driver = webdriver.Chrome(f"{parent_dir}/chromedriver_mac_arm64/chromedriver", options=options)
    driver.get('http://tools.iedb.org/CD4episcore/')

    searchbox = driver.find_element_by_xpath('/html/body/div[3]/form/table/tbody/tr[3]/td[2]/textarea')
    searchbox.send_keys(p)

    sleep(1)

    threshold_button = driver.find_element_by_xpath('/html/body/div[3]/form/table/tbody/tr[9]/td[2]/select/option[10]')
    threshold_button.click()

    sleep(2)

    submit_button = driver.find_elements_by_xpath('/html/body/div[3]/form/table/tbody/tr[12]/th/div/input[1]')
    submit_button[0].click()

    elem = WebDriverWait(driver, 120).until(
        ec.presence_of_element_located((By.XPATH, '/html/body/div[3]/div[1]/h2')))

    for x in range(len(peptide_list)):
        result = driver.find_element_by_xpath(f'/html/body/div[3]/div[1]/div[3]/table/tbody/tr[{x + 1}]/td[6]').text
        e = driver.find_element_by_xpath(f'/html/body/div[3]/div[1]/div[3]/table/tbody/tr[{x + 1}]/td[3]').text
        current_df.loc[current_df["peptide"] == e, "immunogenicity"] = float(result)

    driver.close()


def get_antigenicity(peptide_list):
    options = webdriver.ChromeOptions()
    options.add_argument("--headless=new")
    driver = webdriver.Chrome(f"{parent_dir}/chromedriver_mac_arm64/chromedriver", options=options)
    driver.get('http://www.ddg-pharmfac.net/vaxijen/VaxiJen/VaxiJen.html')

    file_box = driver.find_element_by_xpath('//input[@type="FILE"]')
    file_box.send_keys(peptide_fasta)

    organism_button = driver.find_element_by_xpath('/html/body/div/table/tbody/tr[4]/td[3]/form/table/tbody/tr[2]/td[2]/p/select/option[3]')
    organism_button.click()

    submit_button = driver.find_elements_by_xpath('/html/body/div/table/tbody/tr[4]/td[3]/form/table/tbody/tr[3]/td[2]/input[1]')
    submit_button[0].click()

    sleep(0.5)

    for x in range(len(peptide_list)):
        result = driver.find_element_by_xpath(f'/html/body/div/table/tbody/tr[4]/td[3]/table/tbody/tr/td/b[{3 * (x + 1)}]').text
        e = driver.find_element_by_xpath(f'/html/body/div/table/tbody/tr[4]/td[3]/table/tbody/tr/td/font[{2 * (x + 1)}]').text
        current_df.loc[current_df["peptide"] == e, "antigenicity"] = float(result)

    driver.close()


def get_allergenicity_algpred(peptide_list):
    options = webdriver.ChromeOptions()
    options.add_argument("--headless=new")
    driver = webdriver.Chrome(f"{parent_dir}/chromedriver_mac_arm64/chromedriver", options=options)
    driver.get('https://webs.iiitd.edu.in/raghava/algpred2/batch.html')

    text_box = driver.find_element_by_xpath('/html/body/header/div[3]/section/form/table/tbody/tr/td/font/p/font[1]/textarea')
    text_box.send_keys(pf)

    threshold_value = driver.find_element_by_xpath('/html/body/header/div[3]/section/form/table/tbody/tr/td/font/font/p[3]/font/font[1]/b/select/option[10]')
    threshold_value.click()

    sleep(0.5)

    submit_button = driver.find_elements_by_xpath('/html/body/header/div[3]/section/form/table/tbody/tr/td/font/font/p[3]/font/font[2]/input[2]')
    submit_button[0].click()

    elem = WebDriverWait(driver, 30).until(ec.presence_of_element_located((By.XPATH, '/html/body/header/div[3]/main/h1/strong/font/b')))

    for x in range(len(peptide_list)):
        result = driver.find_element_by_xpath(f'/html/body/header/div[3]/main/div/table[2]/tbody/tr[{x + 1}]/td[5]').text
        current_df.loc[current_df["peptide"] == peptide_list[x], "allergenicity"] = float(result)

    driver.close()


def get_allergenicity_netallergen(peptide_list):
    input_fasta = pf + ">Epitope119" + "\n" + "TIETLMLLALIAAAA"

    options = webdriver.ChromeOptions()
    options.add_argument("--headless=new")
    driver = webdriver.Chrome(f"{parent_dir}/chromedriver_mac_arm64/chromedriver", options=options)
    driver.get('https://services.healthtech.dtu.dk/services/NetAllergen-1.0/')

    try:
        accept_cookies = driver.find_element_by_xpath('/html/body/div[5]/div[4]/div[2]')
    except NoSuchElementException:
        print("No cookies")
    else:
        accept_cookies.click()

    sleep(3)

    text_box = driver.find_element_by_xpath('/html/body/div[3]/main/div/div[3]/div/div[2]/div[1]/form/table/tbody/tr[1]/td/textarea')
    text_box.send_keys(input_fasta)

    sleep(1)

    submit_button = driver.find_element_by_xpath('/html/body/div[3]/main/div/div[3]/div/div[2]/div[1]/form/table/tbody/tr[5]/td/input[1]')
    driver.execute_script("arguments[0].click();", submit_button)

    elem = WebDriverWait(driver, 10000).until(ec.presence_of_element_located((By.XPATH, '/html/body/font/table/tbody/tr/td[3]/h2')))

    output = driver.find_element_by_xpath('/html/body/pre').text
    output = output[:output.find("Explain")]
    table = StringIO(output)
    df = pd.read_table(table, sep=r"\s+")

    count = 0
    for pep in peptide_list:
        result = df["Prediction"][count]
        current_df.loc[current_df["peptide"] == pep, "allergenicity"] = float(result)
        count += 1


def get_protparam(peptide_list):
    for e in peptide_list:
        options = webdriver.ChromeOptions()
        options.add_argument("--headless=new")
        driver = webdriver.Chrome(f"{parent_dir}/chromedriver_mac_arm64/chromedriver", options=options)
        driver.get('https://web.expasy.org/protparam/')

        searchbox = driver.find_element_by_xpath('//*[@id="sib_body"]/form/textarea')
        searchbox.send_keys(e)

        submitButton = driver.find_elements_by_xpath('/html/body/div[2]/div[2]/form/p[1]/input[2]')
        submitButton[0].click()

        results = driver.find_element_by_xpath('/html/body/div[2]/div[2]/pre[2]').text
        pi = results[
             results.find('Theoretical pI: ') + 16:results.find('\n', results.find('Theoretical pI: ') + 16, -1)]
        half_life = results[results.find('The estimated half-life is: ') + 28:results.find('hours', results.find(
            'The estimated half-life is: ') + 28, -1) - 1]
        instability = results[results.find('The instability index (II) is computed to be ') + 45:results.find('\n', results.find('The instability index (II) is computed to be ') + 45, -1)]
        alipathy = results[results.find('Aliphatic index: ') + 17:results.find('\n', results.find('Aliphatic index: ') + 17, -1)]
        gravy = results[results.find('Grand average of hydropathicity (GRAVY): ') + 41:]

        if ">" not in half_life:
            current_df.loc[current_df["peptide"] == e, "half-life"] = float(half_life)
        else:
            # half_life = str(half_life)
            current_df.loc[current_df["peptide"] == e, "half-life"] = float(half_life[half_life.find(">") + 1:])
        current_df.loc[current_df["peptide"] == e, "instability"] = float(instability)
        current_df.loc[current_df["peptide"] == e, "isoelectric point"] = float(pi)
        current_df.loc[current_df["peptide"] == e, "aliphatic index"] = float(alipathy)
        current_df.loc[current_df["peptide"] == e, "GRAVY score"] = float(gravy)

        driver.close()


def get_toxicity(peptide_list):
    options = webdriver.ChromeOptions()
    options.add_argument("--headless=new")
    driver = webdriver.Chrome(f"{parent_dir}/chromedriver_mac_arm64/chromedriver", options=options)
    driver.get('https://webs.iiitd.edu.in/raghava/toxinpred/multi_submit.php')

    submission_box = driver.find_element_by_xpath('//*[@id="input_box"]')
    submission_box.send_keys(pf)

    submit_button = driver.find_elements_by_xpath('/html/body/table[2]/tbody/tr/td/form/fieldset/table[2]/tbody/tr[3]/td/input[2]')
    submit_button[0].click()

    elem = WebDriverWait(driver, 60).until(ec.presence_of_element_located((By.XPATH, '/html/body/div[2]/table/thead/tr[1]/td[1]')))

    for x in range(len(peptide_list)):
        result = driver.find_element_by_xpath(f'/html/body/div[2]/table/tbody/tr[{x + 1}]/td[4]').text
        e = driver.find_element_by_xpath(f'/html/body/div[2]/table/tbody/tr[{x + 1}]/td[2]/a').text
        current_df.loc[current_df["peptide"] == e, "toxicity"] = result

    driver.close()


def get_ifn(peptide_list):
    options = webdriver.ChromeOptions()
    options.add_argument("--headless=new")
    driver = webdriver.Chrome(f"{parent_dir}/chromedriver_mac_arm64/chromedriver", options=options)
    driver.get('https://webs.iiitd.edu.in/raghava/ifnepitope/predict.php')

    submission_box = driver.find_element_by_xpath('//*[@name="sequence"]')
    submission_box.send_keys(pf)

    submit_button = driver.find_elements_by_xpath('//*[@type="submit"]')
    submit_button[0].click()

    elem = WebDriverWait(driver, 120).until(ec.presence_of_element_located((By.XPATH, '//*[@class="sorting_asc"]')))

    entry_number = driver.find_element_by_xpath('/html/body/div/div/div[3]/div/div/div/div/div[1]/label/select/option[4]')
    entry_number.click()

    for x in range(len(peptide_list)):
        result = driver.find_element_by_xpath(f'/html/body/div/div/div[3]/div/div/div/div/table/tbody/tr[{x + 1}]/td[5]').text
        e = driver.find_element_by_xpath(f'/html/body/div/div/div[3]/div/div/div/div/table/tbody/tr[{x + 1}]/td[3]/a').text
        current_df.loc[current_df["peptide"] == e, "IFNg"] = result

    driver.close()


def normalize_data(df_collected_epitopes, mhc):
    immunogenicity_name = "immunogenicity"
    antigenicity_name = "antigenicity"
    allergenicity_name = "allergenicity"
    binding_affinity_name = "binding affinity (nM)"
    immunogenicity_mean = df_collected_epitopes[immunogenicity_name].mean()
    immunogenicity_std = df_collected_epitopes[immunogenicity_name].std()
    antigenicity_mean = df_collected_epitopes[antigenicity_name].mean()
    antigenicity_std = df_collected_epitopes[antigenicity_name].std()
    allergenicity_mean = df_collected_epitopes[allergenicity_name].mean()
    allergenicity_std = df_collected_epitopes[allergenicity_name].std()
    ba_mean = df_collected_epitopes[binding_affinity_name].mean()
    ba_std = df_collected_epitopes[binding_affinity_name].std()
    for i in range(df_collected_epitopes.shape[0]):
        immunogenicity_value = df_collected_epitopes[immunogenicity_name][i]
        antigenicity_value = df_collected_epitopes[antigenicity_name][i]
        allergenicity_value = df_collected_epitopes[allergenicity_name][i]
        ba_value = df_collected_epitopes[binding_affinity_name][i]
        new_immunogenicity_value = (immunogenicity_value - immunogenicity_mean) / immunogenicity_std
        new_antigenicity_value = (antigenicity_value - antigenicity_mean) / antigenicity_std
        new_allergenicity_value = (allergenicity_value - allergenicity_mean) / allergenicity_std
        new_ba_value = (ba_value - ba_mean) / ba_std
        df_collected_epitopes.at[i, "norm immunogenicity"] = new_immunogenicity_value
        df_collected_epitopes.at[i, "norm antigenicity"] = new_antigenicity_value
        df_collected_epitopes.at[i, "norm allergenicity"] = new_allergenicity_value
        df_collected_epitopes.at[i, "norm binding affinity"] = new_ba_value

    immunogenicity_name = "norm immunogenicity"
    antigenicity_name = "norm antigenicity"
    allergenicity_name = "norm allergenicity"
    binding_affinity_name = "norm binding affinity"
    immunogenicity_max = df_collected_epitopes[immunogenicity_name].max()
    immunogenicity_min = df_collected_epitopes[immunogenicity_name].min()
    antigenicity_max = df_collected_epitopes[antigenicity_name].max()
    antigenicity_min = df_collected_epitopes[antigenicity_name].min()
    allergenicity_max = df_collected_epitopes[allergenicity_name].max()
    allergenicity_min = df_collected_epitopes[allergenicity_name].min()
    ba_max = df_collected_epitopes[binding_affinity_name].max()
    ba_min = df_collected_epitopes[binding_affinity_name].min()
    for i in range(df_collected_epitopes.shape[0]):
        immunogenicity_value = df_collected_epitopes[immunogenicity_name][i]
        antigenicity_value = df_collected_epitopes[antigenicity_name][i]
        allergenicity_value = df_collected_epitopes[allergenicity_name][i]
        ba_value = df_collected_epitopes[binding_affinity_name][i]
        if mhc == "I":
            new_immunogenicity_value = (immunogenicity_value - immunogenicity_min)/(immunogenicity_max - immunogenicity_min)
        else:
            new_immunogenicity_value = (immunogenicity_value - immunogenicity_max) / (immunogenicity_min - immunogenicity_max)
        new_antigenicity_value = (antigenicity_value - antigenicity_min)/(antigenicity_max - antigenicity_min)
        new_allergenicity_value = (allergenicity_value - allergenicity_max)/(allergenicity_min - allergenicity_max)
        new_ba_value = (ba_value - ba_max) / (ba_min - ba_max)
        df_collected_epitopes.at[i, immunogenicity_name] = new_immunogenicity_value
        df_collected_epitopes.at[i, antigenicity_name] = new_antigenicity_value
        df_collected_epitopes.at[i, allergenicity_name] = new_allergenicity_value
        df_collected_epitopes.at[i, binding_affinity_name] = new_ba_value
    return df_collected_epitopes


def apply_scoring_function(df_normalized_epitopes, mhc):
    if mhc == "I":
        training_df = pd.read_excel("CD8_training_set.xlsx")
    else:
        training_df = pd.read_excel("CD4_training_set.xlsx")
    X = training_df[["norm immunogenicity", "norm antigenicity", "norm allergenicity", "norm binding affinity"]].values
    y = training_df[["Potential"]]
    regr = linear_model.LinearRegression()
    regr.fit(X, y)

    for i in range(df_normalized_epitopes.shape[0]):
        predicted_potential = regr.predict([[df_normalized_epitopes["norm immunogenicity"][i], df_normalized_epitopes["norm antigenicity"][i], df_normalized_epitopes["norm allergenicity"][i], df_normalized_epitopes["norm binding affinity"][i]]])
        df_normalized_epitopes.at[i, "potential"] = predicted_potential

    # immunogenicity = df_normalized_epitopes["norm immunogenicity"][i]
    # antigenicity = df_normalized_epitopes["norm antigenicity"][i]
    # allergenicity = df_normalized_epitopes["norm allergenicity"][i]
    # ba = df_normalized_epitopes["norm binding affinity"][i]
    # if mhc == "I":
    #     immunogenicity_weight = 0.16824783
    #     antigenicity_weight = -0.01742922
    #     allergenicity_weight = 0.22956087
    #     ba_weight = 0.24842511
    #     predicted_potential = (immunogenicity_weight * immunogenicity) + (antigenicity_weight * antigenicity) + (allergenicity_weight * allergenicity) + (ba_weight * ba)
    #     df_normalized_epitopes.at[i, "potential"] = predicted_potential
    # else:
    #     immunogenicity_weight = 0.11814973
    #     antigenicity_weight = 0.4579201
    #     allergenicity_weight = -0.06767124
    #     ba_weight = 0.18264528
    #     predicted_potential = (immunogenicity_weight * immunogenicity) + (antigenicity_weight * antigenicity) + (allergenicity_weight * allergenicity) + (ba_weight * ba)
    #     df_normalized_epitopes.at[i, "potential"] = predicted_potential
    df_normalized_epitopes_ranked = df_normalized_epitopes.sort_values(by=["potential"], ascending=False)
    return df_normalized_epitopes_ranked


def get_filtered_epitopes(df_ranked_epitopes, mhc):
    top_20_df = df_ranked_epitopes.head(20)
    if mhc == "I":
        df_filtered_epitopes = top_20_df.loc[(top_20_df["half-life"] > 1) & (top_20_df["instability"] < 40) & (top_20_df["toxicity"] == "Non-Toxin")]
    else:
        df_filtered_epitopes = top_20_df.loc[(top_20_df["half-life"] > 1) & (top_20_df["instability"] < 40) & (top_20_df["toxicity"] == "Non-Toxin") & (top_20_df["IFNg"] == "POSITIVE")]
    return df_filtered_epitopes.reset_index(drop=True)


def get_optimized_epitopes(df_filtered_epitopes, mhc):
    filtered_epitopes_str = ""
    op_df = pd.DataFrame(columns=df_filtered_epitopes.columns)
    if mhc == "I":
        java_file = "PopCoverageOptimization.java"
        for i in range(df_filtered_epitopes.shape[0]):
            allele = df_filtered_epitopes["allele"][i]
            peptide = df_filtered_epitopes["peptide"][i]
            filtered_epitopes_str = f"{filtered_epitopes_str}{allele}\t{peptide}\n"
        filtered_epitopes_str = filtered_epitopes_str.rstrip("\n")
        mhc_i = subprocess.run(["java", java_file, filtered_epitopes_str], capture_output=True, text=True)
        output = mhc_i.stdout.rstrip("\n")
        peptide_epitopes_list = output.split("\n")
        for combination in peptide_epitopes_list:
            allele_list = combination[combination.find("HLA"):].split(",")
            pep = combination[:combination.find("HLA") - 1]
            for allele in allele_list:
                index = df_filtered_epitopes.index[(df_filtered_epitopes["allele"] == allele) & (df_filtered_epitopes["peptide"] == pep)][0]
                op_df.loc[len(op_df)] = df_filtered_epitopes.iloc[index]
    else:
        java_file = "CD4PopCoverageOptimization.java"
        allowed_alleles = ["HLA-DMA", "HLA-DMB", "HLA-DPA1", "HLA-DPB1", "HLA-DQA1", "HLA-DQB1", "HLA-DRB1"]
        for i in range(df_filtered_epitopes.shape[0]):
            allele = df_filtered_epitopes["allele"][i]
            allele_check = allele[:allele.find("*")]
            if allele_check not in allowed_alleles or "/" in allele:
                pass
            else:
                peptide = df_filtered_epitopes["peptide"][i]
                filtered_epitopes_str = f"{filtered_epitopes_str}{allele}\t{peptide}\n"
        filtered_epitopes_str = filtered_epitopes_str.rstrip("\n")
        mhc_ii = subprocess.run(["java", java_file, filtered_epitopes_str], capture_output=True, text=True)
        output = mhc_ii.stdout.rstrip("\n")
        peptide_epitopes_list = output.split("\n")
        for combination in peptide_epitopes_list:
            allele_list = combination[combination.find("HLA"):].split(",")
            pep = combination[:combination.find("HLA") - 1]
            for allele in allele_list:
                index = df_filtered_epitopes.index[(df_filtered_epitopes["allele"] == allele) & (df_filtered_epitopes["peptide"] == pep)][0]
                op_df.loc[len(op_df)] = df_filtered_epitopes.iloc[index]
    return op_df


def convert_to_text(df):
    peptide_allele_dict = {}
    for i in range(df.shape[0]):
        peptide = df["peptide"][i]
        allele = df["allele"][i]
        if peptide not in peptide_allele_dict.keys():
            peptide_allele_dict[peptide] = ""
            peptide_allele_dict[peptide] = f"{peptide_allele_dict[peptide]}{allele},"
        else:
            peptide_allele_dict[peptide] = f"{peptide_allele_dict[peptide]}{allele},"
    file_lines = []
    for peptide in peptide_allele_dict.keys():
        alleles = peptide_allele_dict[peptide]
        alleles = alleles[:-1]
        file_lines.append(f"{peptide}\t{alleles}")
    return file_lines


def population_coverage_helper(pop_filename, reg, mhc, pep_filename, plot_folder, path):
    completed_run = subprocess.run(
        ["python", pop_filename, "-p", reg, "-c", mhc, "-f", pep_filename, "--plot", plot_folder], capture_output=True,
        text=True)
    output = completed_run.stdout
    output = output[output.find("World"):]
    output_list = output.split("\n")
    output_list = output_list[:19]
    if path == "pop_cov_path":
        output_list.insert(0, "Region\tCoverage\tAverage Hit\tpc90")
    else:
        output_list.insert(0, "Region\tOptimized Coverage\tAverage Hit\tpc90")
    output = "\n".join(output_list)
    table = StringIO(output)
    df = pd.read_table(table, sep="\t")
    return df


def get_population_coverage(df_filtered_epitopes, df_optimized_epitopes, mhc, c):
    if mhc == "I":
        peptide_allele_filename = f"filtered_epitopes_{c}_mhci.txt"
        optimized_peptide_allele_filename = f"optimized_epitopes_{c}_mhci.txt"
        filtered_epitopes_txt = convert_to_text(df_filtered_epitopes)
        with open(peptide_allele_filename, "w") as fo:
            for line in filtered_epitopes_txt:
                fo.write(f"{line}\n")
        optimized_epitopes_txt = convert_to_text(df_optimized_epitopes)
        with open(optimized_peptide_allele_filename, "w") as fo:
            for line in optimized_epitopes_txt:
                fo.write(f"{line}\n")
    else:
        peptide_allele_filename = f"filtered_epitopes_{c}_mhcii.txt"
        optimized_peptide_allele_filename = f"optimized_epitopes_{c}_mhcii.txt"
        filtered_epitopes_txt = convert_to_text(df_filtered_epitopes)
        with open(peptide_allele_filename, "w") as fo:
            for line in filtered_epitopes_txt:
                fo.write(f"{line}\n")
        optimized_epitopes_txt = convert_to_text(df_optimized_epitopes)
        with open(optimized_peptide_allele_filename, "w") as fo:
            for line in optimized_epitopes_txt:
                fo.write(f"{line}\n")
    new_dir = f"{parent_dir}/Population_Coverage_Results"
    os.makedirs(new_dir, exist_ok=True)
    pop_coverage_filename = f"{parent_dir}/population_coverage/calculate_population_coverage.py"
    plot_output_folder = f"{new_dir}/Plots/Regular_Plots_{c}"
    os.makedirs(plot_output_folder, exist_ok=True)
    optimized_plot_output_folder = f"{new_dir}/Plots/Optimized_Plots_{c}"
    os.makedirs(optimized_plot_output_folder, exist_ok=True)
    regions = ["World", "East Asia", "Northeast Asia", "South Asia", "Southeast Asia", "Southwest Asia", "Europe",
               "East Africa", "West Africa", "Central Africa", "North Africa", "South Africa", "West Indies",
               "North America", "Central America", "South America", "Oceania"]
    region_string = ",".join(regions)
    regular_results = population_coverage_helper(pop_coverage_filename, region_string, mhc_class, peptide_allele_filename, plot_output_folder, "pop_cov_path")
    optimized_results = population_coverage_helper(pop_coverage_filename, region_string, mhc_class, optimized_peptide_allele_filename, optimized_plot_output_folder, "op_pop_cov_path")
    return regular_results, optimized_results


start = time.time()
parent_dir = os.getcwd()
print(parent_dir)
gene_target = "PIK3CA"
gene_file = f"{gene_target}.fasta"
print("Getting gene sequence in fasta format...")
get_gene_sequence(gene_target)
print("Obtained gene sequence")
cancer_mutations_dict = {"CRC": ["R38H", "R88Q", "G106V", "C420R", "E453Q", "E542K", "E545K", "R1023Q", "M1043I", "H1047R"],
                         "Meningioma": ["E110K", "I391M", "R108H", "G914R", "N345K", "E453K", "Y165H", "H1047R", "E545K"],
                         "BC": ["E542K", "E542V", "E545K", "Q546E", "Q546R", "H1047L", "H1047R", "N345K", "E726K", "C420R", "G118D", "E453K", "Q546K", "G1049R", "M1043I", "K111E", "E81K", "E545A", "E545G", "N1044K", "S405P"],
                         "Endometrial": ["E542K", "E542Q", "E545K", "E545G", "G1007R", "Y1021H", "Y1021C", "A1035V", "M1043I", "H1047Y", "H1047R", "G1050D", "T1052K", "H1065L"],
                         "Glioblastoma": ["R88E", "E542K", "E545A", "T1025N", "Y1021N", "R88Q", "P298T", "R310C", "T1031G", "V344G", "E453K", "E545K", "Y1021C", "M1043I", "N1044S", "H1047Y", "G1049S"]}
# cancer_mutations_dict = {"CRC": ["R38H", "E452Q", "M1043I"], "Glioblastoma": ["M1043I"]}
# cancer_mutations_dict = {"CRC": ["R38H", "E542K"], "Glioblastoma": ["E542K"]}
mutations_cancer_dict = {}
df_dict = {}
df_ranked_dict = {}
filtered_dict = {}
# optimized_dict = {}
pop_cov_dict = {}
optimized_pop_cov_dict = {}
cancers = [cancer for cancer in cancer_mutations_dict.keys()]
all_mutations = []
for cancer in cancers:
    for mutation in cancer_mutations_dict[cancer]:
        if mutation not in all_mutations:
            all_mutations.append(mutation)
            mutations_cancer_dict[mutation] = [cancer]
        else:
            mutations_cancer_dict[mutation].append(cancer)
            # df_cancers = pd.DataFrame(my_cancers, columns=["Cancers"])
            # mutations_cancer_dict[mutation] = df_cancers
mhc_classes = ["I", "II"]
# mhc_classes = ["I"]
for mhc_class in mhc_classes:
    print("Making mutant fasta proteins...")
    make_mutant_genes(all_mutations, gene_file)
    print("Done making mutant fasta proteins")
    print("Obtaining all epitopes and binding affinities...")
    mutant_gene_all_epitopes_dict = get_epitopes_ba(all_mutations, mhc_class, mutations_cancer_dict)
    print("Done generating possible epitopes and binding affinities")
    print("Filtering out unmutated epitopes...")
    mutant_gene_mut_epitopes_dict = get_mutant_epitopes(all_mutations, mhc_class, mutant_gene_all_epitopes_dict,
                                                        mutations_cancer_dict)
    print("Done filtering out unmutated epitopes")
    get_peptides_path = f"{parent_dir}/Sequences/getpeptides_autoepicollect.py"
    print("Generating possible peptide sequences for all mutations...")
    get_peptides(all_mutations, mutant_gene_mut_epitopes_dict, mhc_class)
    print("Done generating all peptide sequences")
    if mhc_class == "I":
        final_out = f"all_variables_mhci.xlsx"
    else:
        final_out = f"all_variables_mhcii.xlsx"
    with pd.ExcelWriter(final_out, engine='openpyxl') as w:
        for pm in all_mutations:
            print(f"Getting results for {pm}...")
            immunogenicity_file = f"{parent_dir}/immunogenicity/predict_immunogenicity.py"
            if mhc_class == "I":
                # current_df = pd.read_excel(f"mutant_gene_mut_epitopes_mhci_{cancer}.xlsx", sheet_name=pm)
                current_df = mutant_gene_mut_epitopes_dict[pm].copy()
                peptide_fasta = f"{parent_dir}/Sequences/{pm}peptides_{mhc_class}.fasta"
                with open(peptide_fasta, "r") as f:
                    pf = f.read()
                peptide_file = f"{parent_dir}/Sequences/{pm}peptides_{mhc_class}.txt"
                with open(peptide_file, "r") as fp:
                    p = fp.read()
                with open(peptide_file, "r") as fp:
                    peptides = [line.strip() for line in fp]
                hla_file = "MHCI_HLAs.txt"
                with open(hla_file, "r") as fh:
                    hlas = fh.readlines()
                print(f"Obtaining {pm} immunogenicity...")
                get_local_immunogenicity_mhci()
                print(f"Done obtaining {pm} immunogenicity scores")
                print(f"Obtaining {pm} antigenicity...")
                get_antigenicity(peptides)
                print(f"Done obtaining {pm} antigenicity scores")
                print(f"Obtaining {pm} allergenicity...")
                get_allergenicity_algpred(peptides)
                print(f"Done obtaining {pm} allergenicity predictions")
                print(f"Obtaining {pm} toxicity...")
                get_toxicity(peptides)
                print(f"Done obtaining {pm} toxicity predictions")
                print(f"Obtaining {pm} half-life, instability, isoelectric point, aliphatic index, and GRAVY scores...")
                get_protparam(peptides)
                print(f"Done obtaining {pm} half-life, instability, isoelectric point, aliphatic index, and GRAVY scores")
                print(f"Obtaining {pm} IFN-gamma...")
                get_ifn(peptides)
                print(f"Done obtaining {pm} IFN-gamma predictions")
                df_dict[pm] = current_df
                current_df.to_excel(w, header=True, index=False, sheet_name=pm)
                pd.DataFrame(mutations_cancer_dict[pm], columns=["Cancers"]).to_excel(w, header=True, index=False, sheet_name=pm,
                                                          startcol=current_df.shape[1] + 1, startrow=0)
                print(current_df)
                print(f"Done obtaining all MHC I epitopes and clinical variables for {pm} mutation")
            else:
                # current_df = pd.read_excel(f"mutant_gene_mut_epitopes_mhcii_{cancer}.xlsx", sheet_name=pm)
                current_df = mutant_gene_mut_epitopes_dict[pm].copy()
                peptide_fasta = f"{parent_dir}/Sequences/{pm}peptides_{mhc_class}.fasta"
                with open(peptide_fasta, "r") as f:
                    pf = f.read()
                peptide_file = f"{parent_dir}/Sequences/{pm}peptides_{mhc_class}.txt"
                with open(peptide_file, "r") as fp:
                    p = fp.read()
                with open(peptide_file, "r") as fp:
                    peptides = [line.strip() for line in fp]
                hla_file = "MHCII_HLAs.txt"
                with open(hla_file, "r") as fh:
                    hlas = fh.readlines()
                print(f"Obtaining {pm} immunogenicity...")
                get_immunogenicity_mhcii(peptides)
                print(f"Done obtaining {pm} immunogenicity scores")
                print(f"Obtaining {pm} antigenicity...")
                get_antigenicity(peptides)
                print(f"Done obtaining {pm} antigenicity scores")
                print(f"Obtaining {pm} allergenicity...")
                get_allergenicity_netallergen(peptides)
                # get_allergenicity_algpred(peptides)
                print(f"Done obtaining {pm} allergenicity predictions")
                print(f"Obtaining {pm} toxicity...")
                get_toxicity(peptides)
                print(f"Done obtaining {pm} toxicity predictions")
                print(f"Obtaining {pm} half-life, instability, isoelectric point, aliphatic index, and GRAVY scores...")
                get_protparam(peptides)
                print(f"Done obtaining {pm} half-life, instability, isoelectric point, aliphatic index, and GRAVY scores")
                print(f"Obtaining {pm} IFN-gamma...")
                get_ifn(peptides)
                print(f"Done obtaining {pm} IFN-gamma predictions")
                df_dict[pm] = current_df
                current_df.to_excel(w, header=True, index=False, sheet_name=pm)
                pd.DataFrame(mutations_cancer_dict[pm], columns=["Cancers"]).to_excel(w, header=True, index=False, sheet_name=pm,
                                                          startcol=current_df.shape[1] + 1, startrow=0)
                print(current_df)
                print(f"Done obtaining all MHC II epitopes and clinical variables for {pm} mutation")
    print("Done obtaining all epitopes and clinical variables for all mutations")

    if mhc_class == "I":
        ranked_final_out = "all_variables_ranked_mhci.xlsx"
    else:
        ranked_final_out = "all_variables_ranked_mhcii.xlsx"
    with pd.ExcelWriter(ranked_final_out, engine='openpyxl') as w:
        for pm in all_mutations:
            print(f"Normalizing and ranking epitopes for {pm} mutation...")
            current_df = df_dict[pm]
            normalized_df = normalize_data(current_df, mhc_class)
            ranked_df = apply_scoring_function(normalized_df, mhc_class)
            df_ranked_dict[pm] = ranked_df
            ranked_df.to_excel(w, header=True, index=False, sheet_name=pm)
            pd.DataFrame(mutations_cancer_dict[pm], columns=["Cancers"]).to_excel(w, header=True, index=False,
                                                                                  sheet_name=pm,
                                                                                  startcol=current_df.shape[1] + 1,
                                                                                  startrow=0)
            print(f"Done normalizing and ranking epitopes for {pm} mutation")
        print("Done normalizing and ranking epitopes for all mutations")

    if mhc_class == "I":
        top_epitopes_out = "top_epitopes_mhci.xlsx"
    else:
        top_epitopes_out = "top_epitopes_mhcii.xlsx"
    with pd.ExcelWriter(top_epitopes_out, engine='openpyxl') as w:
        for pm in all_mutations:
            print(f"Filtering epitopes for {pm} mutation...")
            ranked_df = df_ranked_dict[pm]
            filtered_df = get_filtered_epitopes(ranked_df, mhc_class)
            filtered_dict[pm] = filtered_df
            filtered_df.to_excel(w, header=True, index=False, sheet_name=pm)
            pd.DataFrame(mutations_cancer_dict[pm], columns=["Cancers"]).to_excel(w, header=True, index=False,
                                                                                  sheet_name=pm,
                                                                                  startcol=filtered_df.shape[1] + 1,
                                                                                  startrow=0)
            print(f"Done filtering epitopes for {pm} mutation")
    if mhc_class == "I":
        top_filtered_epitopes_out = "top_epitopes_by_cancer_mhci.xlsx"
    else:
        top_filtered_epitopes_out = "top_epitopes_by_cancer_mhcii.xlsx"
    with pd.ExcelWriter(top_filtered_epitopes_out, engine='openpyxl') as w:
        for cancer in cancers:
            sum_filtered_epitopes_df = pd.DataFrame()
            for pm in cancer_mutations_dict[cancer]:
                filtered_df = filtered_dict[pm]
                sum_filtered_epitopes_df = pd.concat([sum_filtered_epitopes_df, filtered_df], ignore_index=True)
            pop_cov_dict[cancer] = sum_filtered_epitopes_df
            sum_filtered_epitopes_df.to_excel(w, header=True, index=False, sheet_name=cancer)
    print("Done filtering epitopes for all mutations")

    print(f"Obtaining MHC Class {mhc_class} population coverage results...")
    # if mhc_class == "I":
    #     op_top_epitopes_out = "optimized_top_epitopes_mhci.xlsx"
    # else:
    #     op_top_epitopes_out = "optimized_top_epitopes_mhcii.xlsx"
    # with pd.ExcelWriter(op_top_epitopes_out, engine='openpyxl') as w:
    #     for pm in all_mutations:
    #         print(f"Optimizing epitopes for {pm} mutation...")
    #         filtered_df = filtered_dict[pm]
    #         optimized_df = get_optimized_epitopes(filtered_df, mhc_class)
    #         optimized_dict[pm] = optimized_df
    #         optimized_df.to_excel(w, header=True, index=False, sheet_name=pm)
    #         pd.DataFrame(mutations_cancer_dict[pm], columns=["Cancers"]).to_excel(w, header=True, index=False,
    #                                                                               sheet_name=pm,
    #                                                                               startcol=optimized_df.shape[1] + 1,
    #                                                                               startrow=0)
    #         print(f"Done optimizing epitopes for {pm} mutation")
    if mhc_class == "I":
        top_op_epitopes_out = "optimized_epitopes_by_cancer_mhci.xlsx"
    else:
        top_op_epitopes_out = "optimized_epitopes_by_cancer_mhcii.xlsx"
    with pd.ExcelWriter(top_op_epitopes_out, engine='openpyxl') as w:
        for cancer in cancers:
            print(f"Optimizing epitopes for {cancer}...")
            sum_filtered_epitopes_df = pop_cov_dict[cancer]
            if sum_filtered_epitopes_df.empty:
                print(f"No filtered MHC Class {mhc_class} epitopes for {cancer}, cannot perform population coverage analysis")
            else:
                optimized_df = get_optimized_epitopes(sum_filtered_epitopes_df, mhc_class)
                optimized_pop_cov_dict[cancer] = optimized_df
                optimized_df.to_excel(w, header=True, index=False, sheet_name=cancer)
                print(f"Done optimizing epitopes for {cancer}")
    print("Done optimizing epitopes for all cancers")

    if mhc_class == "I":
        population_coverage_out = "population_results_by_cancer_mhci.xlsx"
    else:
        population_coverage_out = "population_results_by_cancer_mhcii.xlsx"
    with pd.ExcelWriter(population_coverage_out, engine='openpyxl') as w:
        for cancer in cancers:
            sum_filtered_epitopes_df = pop_cov_dict[cancer]
            if not sum_filtered_epitopes_df.empty:
                print(f"Obtaining population coverage for {cancer}...")
                optimized_df = optimized_pop_cov_dict[cancer]
                regular_population_coverage, optimized_population_coverage = get_population_coverage(sum_filtered_epitopes_df, optimized_df, mhc_class, cancer)
                regular_population_coverage.to_excel(w, header=True, index=False, sheet_name=cancer)
                optimized_population_coverage.to_excel(w, header=True, index=False, sheet_name=cancer, startcol=regular_population_coverage.shape[1] + 1, startrow=0)
                print(f"Done obtaining population coverage for {cancer}")
    print(f"Done obtaining MHC Class {mhc_class} population coverage results for all cancers")
print("AutoEpiCollect complete, please click the button to the right to see your epitopes for each mutation")
print("All other outputted Excel spreadsheets will be in the same directory as this program")






