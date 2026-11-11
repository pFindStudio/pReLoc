import pandas as pd
import numpy as np
import itertools
import copy
from generate_decoy import generate_decoy
from process_pfindres import process_pfindres
from process_pChemres import get_deltamass
from process_MGF import process_MGF
import sys
import os


def modification_parser(modification_ini_path):
    modification_dict = {}
    with open(modification_ini_path, "r", encoding="utf-8") as f:
        modification_lines = f.readlines()  # Read the complete modification.ini file.

    for index, line in enumerate(modification_lines):
        if "@NUMBER" in line or line == "\n":  # @NUMBER_MODIFICATION=1610
            continue
        if line[0:4] == "name":  # name1=2-dimethylsuccinyl[C] 1
            # 2-dimethylsuccinyl[C]=C NORMAL 144.042259 144.042259 0 H(8)C(6)O(4)
            t_line = modification_lines[index + 1]
            mod_name = t_line.split("=")[0]  # 2-dimethylsuccinyl[C]

            mod_site = t_line.split("=")[1].split(" ")[0]  # C
            mod_mass = float(t_line.split("=")[1].split(" ")[2])  # 144.042259
            if t_line.split("=")[1].split(" ")[4] == "0":
                mod_loss = 0  # mod_loss=0
            else:
                mod_loss = float(t_line.split("=")[1].split(" ")[5])
            modification_dict[mod_name] = {}  # {'2-dimethylsuccinyl[C]': {}}
            if (
                "AnyN-term" in mod_name[mod_name.find("[") + 1 : mod_name.find("]")]
            ) or (
                "ProteinN-term" in mod_name[mod_name.find("[") + 1 : mod_name.find("]")]
            ):
                modification_dict[mod_name]["mod_site_restrict"] = 0
            elif (
                "AnyC-term" in mod_name[mod_name.find("[") + 1 : mod_name.find("]")]
            ) or (
                "ProteinC-term" in mod_name[mod_name.find("[") + 1 : mod_name.find("]")]
            ):
                modification_dict[mod_name]["mod_site_restrict"] = -1
            else:
                modification_dict[mod_name]["mod_site_restrict"] = 1

            modification_dict[mod_name]["mod_site"] = mod_site
            modification_dict[mod_name]["mod_mass"] = mod_mass
            modification_dict[mod_name]["mod_loss"] = mod_loss
    return modification_dict


def get_specname_list(param):
    specname_list = []
    with open(param.MGFPath, "r") as MGF:
        for line in MGF:
            if line.startswith("TITLE"):
                line = line.strip()
                specname_list.append(line.split("=")[1])
    return specname_list


def process_PEAKSRes(param, modification_dict):
    # modification_dict = modification_parser(modification_ini_path)
    Mod_dict = {}
    for modname in param.ModName:
        modmass = modification_dict[modname]["mod_mass"]
        if int(modmass * 100) / 100.0 != round(modmass, 2):
            if (
                int(modmass * 100) / 100.0 in Mod_dict.keys()
                or round(modmass, 2) in Mod_dict.keys()
            ):
                KeyError("Peaks Mod Mass Repeat!")
            Mod_dict[int(modmass * 100) / 100.0] = modname
            Mod_dict[round(modmass, 2)] = modname
        else:
            if int(modmass * 100) / 100.0 in Mod_dict.keys():
                KeyError("Peaks Mod Mass Repeat!")
            Mod_dict[int(modmass * 100) / 100.0] = modname

    specname_list = get_specname_list(param)

    peaksRes_df = pd.read_csv(param.InputPath)
    with open(param.AutoCustomInputPath, "w") as oc:
        oc.write("spec_name\tpep\tmod\tscore\n")
        for index, row in peaksRes_df.iterrows():
            pep = row["Peptide"]
            spec_name = specname_list[row["Scan"] - 1]
            score = row["ALC (%)"]

            mod_index = 0
            i = 0
            new_pep = ""
            mod = ""
            while i < len(pep):
                if pep[i] == "(":
                    mass = float(pep[i + 2 : i + 7])

                    if mass in Mod_dict.keys():
                        if "N-term" in Mod_dict[mass]:
                            mod += str(0) + "," + Mod_dict[mass] + ";"
                        elif "C-term" in Mod_dict[mass]:
                            mod += str(mod_index + 1) + "," + Mod_dict[mass] + ";"
                        else:
                            mod += str(mod_index) + "," + Mod_dict[mass] + ";"
                    else:
                        KeyError("Mod mass error!")
                    i = i + 8
                    continue
                if pep[i] == "I":
                    new_pep += "L"
                else:
                    new_pep += pep[i]
                mod_index += 1
                i += 1
            oc.write(spec_name + "\t" + new_pep + "\t" + mod + "\t" + str(score) + "\n")
    return


def process_pFindspectra(param):
    pFindspectra_df = pd.read_csv(param.InputPath, sep="\t")
    with open(param.AutoCustomInputPath, "w") as oc:
        oc.write("spec_name\tpep\tmod\tscore\n")
        for index, row in pFindspectra_df.iterrows():
            spec_name = row["File_Name"]
            pep = row["Sequence"].replace("I", "L")
            mod = row["Modification"]
            if isinstance(mod, float) and np.isnan(mod):
                mod = ""
            score = -np.log10(float(row["Final_Score"]))
            oc.write(spec_name + "\t" + pep + "\t" + mod + "\t" + str(score) + "\n")
    return


def process_modpep(modpep):
    leftc = False
    iner_leftc = False

    pep = ""
    mod = ""
    modname = ""
    modsite = ""
    for i in range(len(modpep)):
        if modpep[i] == "_":
            continue

        if not leftc and not iner_leftc and modpep[i] == "(":
            leftc = True
            continue
        if leftc and not iner_leftc and modpep[i] == "(":
            iner_leftc = True
            continue
        if leftc and iner_leftc and modpep[i] == ")":
            iner_leftc = False
            continue
        if leftc and not iner_leftc and modpep[i] == ")":
            leftc = False
            if "N-term" in modsite:
                mod += "0," + modname + "[" + modsite + "];"
            elif "C-term" in modsite:
                mod += str(len(pep) + 1) + "," + modname + "[" + modsite + "];"
            else:
                mod += str(len(pep)) + "," + modname + "[" + pep[-1] + "];"
            modname = ""
            modsite = ""
            continue

        if not leftc and not iner_leftc:
            if modpep[i] == "I":
                pep += "L"
            else:
                pep += modpep[i]
            continue
        if leftc and not iner_leftc and modpep[i] != " ":
            modname += modpep[i]
            continue
        if leftc and iner_leftc and modpep[i] != " ":
            modsite += modpep[i]
            continue
    return pep, mod


def process_MQRes(param):
    MQRes_df = pd.read_csv(param.InputPath, sep="\t")
    with open(param.AutoCustomInputPath, "w") as oc:
        oc.write("spec_name\tpep\tmod\tscore\n")
        for index, row in MQRes_df.iterrows():
            spec_name = (
                str(row["Raw file"])
                + "."
                + str(row["Scan number"])
                + "."
                + str(row["Scan number"])
                + "."
                + str(row["Charge"])
                + "."
                + str(row["Scan event number"])
                + ".dta"
            )
            modpep = row["Modified sequence"]
            pep, mod = process_modpep(modpep)
            score = row["Score"]
            oc.write(spec_name + "\t" + pep + "\t" + mod + "\t" + str(score) + "\n")
    return


def combine_unique_tuples(dictionary):
    result = []
    keys = list(dictionary.keys())
    for i in range(len(keys)):
        for j in range(i + 1, len(keys)):
            for tuple1 in dictionary[keys[i]]:
                for tuple2 in dictionary[keys[j]]:
                    if isinstance(tuple1, tuple) and isinstance(tuple2, tuple):
                        if len(set(tuple1).intersection(set(tuple2))) == 0:
                            result.append((keys[i], tuple1) + (keys[j], tuple2))
    return result


def reLocate_CInputFile(param):
    scoremodaa = {}
    # scoremodnum = {}
    for scoreMod in param.scoreModName:
        mt = scoreMod.split("[")[0]
        ms = scoreMod.split("[")[1].split("]")[0].replace("I", "L")
        postfix = "]".join(scoreMod.split("\t")[0].split("]")[1:])
        if mt not in scoremodaa.keys():
            scoremodaa[mt] = {}
            scoremodaa[mt]["aa"] = []
            scoremodaa[mt]["postfix"] = []
            # scoremodnum[mt] = {}
        if ms not in scoremodaa[mt]:
            scoremodaa[mt]["aa"].append(ms)
            scoremodaa[mt]["postfix"].append(postfix)
            # scoremodnum[mt][aa] = 0

    if not param.ispFind:
        AutoCus_df = pd.read_csv(param.InputPath, sep="\t")
    else:
        AutoCus_df = pd.read_csv(param.AutoCustomInputPath, sep="\t")
    has_relocation_candidates = False
    with open(param.AutoCustomInputPath, "w") as oc, open(
        param.DirectResPath, "w"
    ) as dc:
        oc.write("spec_name\tseq\tmod\tscore\n")
        dc.write("spec_name\tseq\tori_mod\tmod\tscoremodsite_num\n")
        for index, row in AutoCus_df.iterrows():
            spec_name = row["spec_name"]
            pep = row["seq"].replace("I", "L")
            mod = row["mod"]
            score = row["score"]

            modnum = {}
            modtype = []
            if isinstance(mod, float) and np.isnan(mod):
                continue
            except_lst = []
            for modinfo in mod.split(";")[:-1]:
                mt = modinfo.split(",")[1].split("[")[0]
                ma = modinfo.split(",")[1].split("[")[1].split("]")[0].replace("I", "L")
                if mt not in modnum.keys():
                    modnum[mt] = 0
                    modtype.append(mt)
                if mt not in scoremodaa.keys():
                    mindex = int(modinfo.split(",")[0])-1
                    if mindex >= 0 and mindex < len(pep):
                        except_lst.append(mindex)
                modnum[mt] += 1

            combs_dict = {}
            mt_list = []
            for mt in modtype:
                modsites = []
                if mt in scoremodaa.keys():
                    mt_list.append(mt)
                    for x in scoremodaa[mt]["aa"]:
                        if "N-term" in x:
                            modsites.append(0)
                            break
                    for i in range(len(pep)):
                        if pep[i] in scoremodaa[mt]["aa"] and i not in except_lst:
                            modsites.append(i)
                    for x in scoremodaa[mt]["aa"]:
                        if "C-term" in x:
                            modsites.append(len(pep) - 1)
                            break
                    combs_dict[mt] = list(itertools.combinations(modsites, modnum[mt]))
                    

            combs = []
            if len(combs_dict) > 1:
                combs = combine_unique_tuples(combs_dict)
            else:
                combs = [
                    (key, value) for key in combs_dict for value in combs_dict[key]
                ]

            modstring_ = []
            flag = True
            for modinfo in mod.split(";")[:-1]:
                for mt in mt_list:
                    if mt in modinfo:
                        flag = False
                        break
                    else:
                        flag = True
                if flag:
                    modstring_.append(modinfo + ";")
            # if len(combs) == 1:
            for comb in combs:
                modstring = copy.deepcopy(modstring_)
                for i in range(int(len(comb) / 2)):
                    mt = comb[i * 2]
                    sites = comb[i * 2 + 1]
                    for site in sites:
                        if pep[site] in scoremodaa[mt]["aa"]:
                            postindex = scoremodaa[mt]["aa"].index(pep[site])
                            modstring.append(
                                str(site + 1)
                                + ","
                                + mt
                                + "["
                                + pep[site]
                                + "]"
                                + scoremodaa[mt]["postfix"][postindex]
                                + ";"
                            )
                        elif site == 0:
                            for aa_i in scoremodaa[mt]["aa"]:
                                if "N-term" in aa_i:
                                    postindex = scoremodaa[mt]["aa"].index(aa_i)
                                    modstring.append(
                                        str(0)
                                        + ","
                                        + mt
                                        + "["
                                        + aa_i
                                        + "]"
                                        + scoremodaa[mt]["postfix"][postindex]
                                        + ";"
                                    )
                                    break
                        elif site == len(pep) - 1:
                            for aa_i in scoremodaa[mt]["aa"]:
                                if "C-term" in aa_i:
                                    postindex = scoremodaa[mt]["aa"].index(aa_i)
                                    modstring.append(
                                        str(len(pep))
                                        + ","
                                        + mt
                                        + "["
                                        + aa_i
                                        + "]"
                                        + scoremodaa[mt]["postfix"][postindex]
                                        + ";"
                                    )
                                    break
                modstring = sorted(modstring, key=lambda x: int(x.split(",")[0]))
                # Concatenate the annotations into one string.
                new_mod = "".join(modstring)
                if len(combs) > 1:
                    has_relocation_candidates = True
                if len(combs) == 1:
                    scoremodsite_num = sum(modnum[mt] for mt in mt_list)
                    dc.write(
                        spec_name
                        + "\t"
                        + pep
                        + "\t"
                        + mod
                        + "\t"
                        + new_mod
                        + "\t"
                        + str(scoremodsite_num)
                        + "\n"
                    )
                elif new_mod == mod:
                    oc.write(
                        spec_name
                        + "\t"
                        + pep
                        + "\t"
                        + new_mod
                        + "\t"
                        + str(score)
                        + "\n"
                    )
                else:
                    oc.write(
                        spec_name + "\t" + pep + "\t" + new_mod + "\t" + str(-1) + "\n"
                    )
    return has_relocation_candidates

def resource_path(relative_path):
    try:
        # PyInstaller stores the bundled resource directory in _MEIPASS.
        base_path = sys._MEIPASS
    except Exception:
        base_path = os.path.dirname(os.path.abspath(__file__))
    
    return os.path.join(base_path, relative_path)

def init_processOtherRes(param):
    scoremodnum = {}
    modification_dict = modification_parser(resource_path(os.path.join("local_denovo", "modification.ini")))
    for x in param.ModName:
        if "\t" in x:
            modname = x.split("\t")[0]
            modmass = x.split("\t")[1]
            modaa = modname.split("[")[1].split("]")[0].replace("I", "L")
            modification_dict[modname] = {}
            modification_dict[modname]["mod_site"] = modaa
            modification_dict[modname]["mod_mass"] = modmass
            modification_dict[modname]["mod_loss"] = 0
    if param.ispFind and not param.ispChem:
        process_pfindres(param)
    if param.ispChem:
        get_deltamass(param)
    process_MGF(param)

    has_relocation_candidates = True
    if param.score_isMod:
        has_relocation_candidates = reLocate_CInputFile(param)
    if param.isFLRest and has_relocation_candidates:
        generate_decoy(param, modification_dict)
    return has_relocation_candidates
