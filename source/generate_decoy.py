import pandas as pd
import random
from write_param import AUTO_MOD_WARNING


DECOY_RANDOM_SEED = 0

LEGACY_AUTO_MOD_WARNING = "# ispfindOpen=true \u65f6\uff0c\u672c\u533a\u6bb5\u4f1a\u81ea\u52a8\u8986\u76d6\uff0c\u8bf7\u52ff\u624b\u52a8\u5199\u5165\u4fee\u9970\u6761\u76ee"


def not_minus_one_to_one(x):
    return 1 if x != -1 else 0


def get_targetdf(target10):
    target10_df = pd.read_csv(target10, sep="\t")
    # target10_df.rename(columns={"score": "target_label"}, inplace=True)
    target10_df["score"] = target10_df["score"].apply(not_minus_one_to_one)
    return target10_df


def swap_char_at_index(s, index, positions_to_exclude):
    positions_to_exclude = set(positions_to_exclude)
    index = index - 1
    # if index < 0 or index >= len(s) or index in positions_to_exclude:
    #     return s, index

    available_positions = [
        i for i in range(len(s)) if i != index and i not in positions_to_exclude
    ]

    if not available_positions:
        return s, index

    swap_index = random.choice(available_positions)
    s_list = list(s)
    s_list[index], s_list[swap_index] = s_list[swap_index], s_list[index]
    swapped_s = "".join(s_list)

    return swapped_s, (swap_index + 1)


def update_and_sort_modinfo(modinfo_str, old_pos, new_pos):
    modinfo_list = modinfo_str.split(";")[:-1]

    updated_modinfo_list = []
    for modinfo in modinfo_list:
        pos, mod = modinfo.split(",")
        if int(pos) == old_pos:
            pos = str(new_pos)
        elif int(pos) == new_pos:
            pos = str(old_pos)
        updated_modinfo_list.append((int(pos), mod))

    sorted_modinfo_list = sorted(updated_modinfo_list, key=lambda x: x[0])

    sorted_modinfo_str = (
        ";".join([f"{pos},{mod}" for pos, mod in sorted_modinfo_list]) + ";"
    )
    return sorted_modinfo_str


def get_decoy1(target10_df, scoremodname_list, scoremodaa_dict):
    new_df = target10_df.copy()
    decoy_list = []
    for index, row in target10_df.iterrows():
        spec_name = row["spec_name"]
        pep = row["seq"]
        mod = row["mod"]

        phomod_list = []
        othermod_list = []
        modinfo_list = mod.split(";")[:-1]
        for modinfo in modinfo_list:
            modindex = int(modinfo.split(",")[0])
            modname = modinfo.split(",")[1]
            if modname.split("[")[0] in scoremodname_list:
                phomod_list.append((modindex, modname.split("[")[0]))
            else:
                othermod_list.append((modindex, modname.split("[")[0]))

        positions_by_mod = {}
        for key, value in scoremodaa_dict.items():
            lst = []
            for i in range(len(pep)):
                if pep[i] in value:
                    lst.append(i)
            positions_by_mod[key] = lst

        for phomod in phomod_list:
            new_pep, new_phoindex = swap_char_at_index(
                pep, phomod[0], positions_by_mod[phomod[1]]
            )
            new_mod = update_and_sort_modinfo(mod, phomod[0], new_phoindex)
            new_decoy = {
                "spec_name": spec_name,
                "seq": new_pep,
                "mod": new_mod,
                "score": -1,
            }
            decoy_list.append(new_decoy)
    decoy_df = pd.DataFrame(decoy_list)
    new_df = pd.concat([new_df, decoy_df], ignore_index=True)
    return new_df


def get_new_mod(modinfo_str, pepindex_lst, scoremodaa_lst, ori_index, pep,scoremodname_list):
    lst_ = [item for item in pepindex_lst if item not in scoremodaa_lst]
    # random_index = random.choice(lst)

    thresh = 4
    lst = []
    if lst_[-1] + 1 - ori_index > thresh or lst_[0] + 1 - ori_index < -thresh:
        for item in pepindex_lst:
            if item + 1 - ori_index > thresh or item + 1 - ori_index < -thresh:
                lst.append(item)
    else:
        lst = lst_
    random_index = random.choice(lst)
    # if lst[-1] + 1 - ori_index > 3 or ori_index - lst[0] - 1 > 3:
    #     while random_index + 1 - ori_index <= 3 and random_index + 1 - ori_index >= -3:
    #         random_index = random.choice(lst)

    modinfo_list = modinfo_str.split(";")[:-1]

    updated_modinfo_list = []
    for modinfo in modinfo_list:
        pos, mod = modinfo.split(",")
        if int(pos) == ori_index and mod.split("[")[0] in scoremodname_list:
            pos = str(random_index + 1)
            modname = mod.split("[")[0]
            mod = modname + "[" + pep[random_index] + "]"
        elif int(pos) == random_index + 1:
            KeyError("Deocy Generate Error!")
        updated_modinfo_list.append((int(pos), mod))

    sorted_modinfo_list = sorted(updated_modinfo_list, key=lambda x: x[0])

    sorted_modinfo_str = (
        ";".join([f"{pos},{mod}" for pos, mod in sorted_modinfo_list]) + ";"
    )
    return sorted_modinfo_str


def get_decoy2(target10_df, scoremodname_list, scoremodaa_dict):
    new_df = target10_df.copy()
    deocy_list = []
    for index, row in target10_df.iterrows():
        spec_name = row["spec_name"]
        pep = row["seq"]
        mod = row["mod"]

        phomod_list = []
        othermod_list = []
        othermodindex_list = []
        modinfo_list = mod.split(";")[:-1]
        for modinfo in modinfo_list:
            modindex = int(modinfo.split(",")[0])
            modname = modinfo.split(",")[1]
            if modname.split("[")[0] in scoremodname_list:
                phomod_list.append((modindex, modname.split("[")[0]))
            else:
                othermod_list.append((modindex, modname.split("[")[0]))
                othermodindex_list.append(modindex - 1)

        scoremodaa_dict1 = {}
        for key, value in scoremodaa_dict.items():
            lst = []
            for i in range(len(pep)):
                if pep[i] in value:
                    lst.append(i)
            scoremodaa_dict1[key] = lst + othermodindex_list

        pepindex_lst = list(range(len(pep)))
        for phomod in phomod_list:
            new_mod = get_new_mod(
                mod, pepindex_lst, scoremodaa_dict1[phomod[1]], phomod[0], pep,scoremodname_list
            )
            new_decoy = {
                "spec_name": spec_name,
                "seq": pep,
                "mod": new_mod,
                "score": -1,
            }
            deocy_list.append(new_decoy)
    decoy_df = pd.DataFrame(deocy_list)
    new_df = pd.concat([new_df, decoy_df], ignore_index=True)
    return new_df


def write_new_param(
    paramFilePath, CparamFilePath, modname_list, scoremodaa_dict, modification_dict
):
    writeline = ""
    line_dict = {}
    linename = ""
    aa_lst = [
        "A",
        "R",
        "N",
        "D",
        "C",
        "E",
        "Q",
        "G",
        "H",
        "I",
        "L",
        "K",
        "M",
        "F",
        "P",
        "S",
        "T",
        "W",
        "Y",
        "V",
    ]
    with open(paramFilePath, "r", encoding="utf-8-sig") as pF:
        for line in pF.readlines():
            if line.strip() == LEGACY_AUTO_MOD_WARNING:
                line = AUTO_MOD_WARNING + "\n"
            flag = True
            if line.startswith("["):
                linename = line[1:-2]
            for x in modname_list:
                if x in line and (linename == "mods" or linename == "scoreModName"):
                    flag = False
            if flag and line != "\n":
                if linename in line_dict.keys():
                    line_dict[linename] += line
                elif linename not in line_dict.keys():
                    line_dict[linename] = line
    with open(CparamFilePath, "w", encoding="utf-8") as CpF:
        for key, value in line_dict.items():
            if key == "mods":
                writeline += value
                mod_indices = []
                for line in value.splitlines():
                    entry_name = line.split("=", 1)[0]
                    if entry_name.startswith("mod") and entry_name[3:].isdigit():
                        mod_indices.append(int(entry_name[3:]))
                index = max(mod_indices, default=0)
                indexname = "mod"
                for x in modname_list:
                    for aa in aa_lst:
                        index = index + 1
                        modname_all = x + "[" + scoremodaa_dict[x][0] + "]"
                        writeline += (
                            indexname
                            + str(index)
                            + "="
                            + x
                            + "["
                            + aa
                            + "]\t"
                            + str(modification_dict[modname_all]["mod_mass"])
                            + "\n"
                        )
                writeline += "\n"
            elif key == "scoreModName":
                writeline += value
                index = 0
                indexname = "ModName"
                for x in modname_list:
                    for aa in aa_lst:
                        index = index + 1
                        modname_all = x + "[" + scoremodaa_dict[x][0] + "]"
                        writeline += (
                            indexname + str(index) + "=" + x + "[" + aa + "]" + "\n"
                        )
                writeline += "\n"
            else:
                writeline += value
                writeline += "\n"
        CpF.write(writeline)
    return


def generate_decoy(param, modification_dict):
    random.seed(DECOY_RANDOM_SEED)
    target10_df = get_targetdf(param.AutoCustomInputPath)
    scoremodname_list = []
    scoremodaa_dict = {}
    for scoremod in param.scoreModName:
        modname = scoremod.split("[")[0]
        modaa = scoremod.split("[")[1].split("]")[0]
        if modname not in scoremodname_list:
            scoremodname_list.append(modname)
            scoremodaa_dict[modname] = [modaa]
        else:
            scoremodaa_dict[modname].append(modaa)
    if param.isFLRest == 1:
        target_decoy_df = get_decoy1(target10_df, scoremodname_list, scoremodaa_dict)
    elif param.isFLRest == 2:
        target_decoy_df = get_decoy2(target10_df, scoremodname_list, scoremodaa_dict)
        if not param.ispChem:
            write_new_param(
                param.paramFilePath,
                param.CparamFilePath,
                scoremodname_list,
                scoremodaa_dict,
                modification_dict,
            )
    target_decoy_df = (
        target_decoy_df[["spec_name", "seq", "mod", "score"]]
        .drop_duplicates()
        .reset_index(drop=True)
    )
    target_decoy_df.to_csv(param.AutoCustomInputPath, sep="\t", index=False)
