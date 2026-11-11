from resources import resource_path
import csv


AUTO_MOD_WARNING = (
    "# When ispfindOpen=true, this section is overwritten automatically; "
    "do not add modification entries manually"
)


def _replace_mod_entries(lines, section_name, prefix, modifications):
    start = next(
        (i for i, line in enumerate(lines) if line.strip() == section_name), None
    )
    if start is None:
        raise ValueError(f"Parameter file is missing the {section_name} section")
    end = next(
        (
            i
            for i in range(start + 1, len(lines))
            if lines[i].strip().startswith("[")
        ),
        len(lines),
    )

    block = lines[start + 1 : end]
    if AUTO_MOD_WARNING not in [line.strip() for line in block]:
        block.insert(0, AUTO_MOD_WARNING + "\n")

    entry_indices = [
        i
        for i, line in enumerate(block)
        if "=" in line and not line.lstrip().startswith("#")
    ]
    entries = [
        f"{prefix}{index}={modification}\n"
        for index, modification in enumerate(modifications, 1)
    ]

    if entry_indices:
        insert_at = entry_indices[0]
        entry_index_set = set(entry_indices)
        block = [
            line for i, line in enumerate(block) if i not in entry_index_set
        ]
        block[insert_at:insert_at] = entries
    else:
        trailing_blank_lines = []
        while block and not block[-1].strip():
            trailing_blank_lines.insert(0, block.pop())
        block.extend(entries)
        block.extend(trailing_blank_lines)

    return lines[: start + 1] + block + lines[end:]


def update_pfind_modifications(param_path, pfind_path):
    modifications = {"[mods]": [], "[NTerm_mods]": [], "[CTerm_mods]": []}
    seen = set()

    with open(pfind_path, "r", encoding="utf-8-sig", newline="") as pfind_file:
        reader = csv.DictReader(pfind_file, delimiter="\t")
        if not reader.fieldnames or "Modification" not in reader.fieldnames:
            raise ValueError("The pFind result file is missing the Modification column")
        for row in reader:
            for item in (row["Modification"] or "").split(";"):
                item = item.strip()
                if not item:
                    continue
                if "," not in item or "[" not in item or "]" not in item:
                    raise ValueError(f"Cannot parse Modification entry: {item}")
                mod_name = item.split(",", 1)[1].strip()
                if mod_name in seen:
                    continue
                seen.add(mod_name)
                mod_site = mod_name.split("[", 1)[1].split("]", 1)[0]
                if "N-term" in mod_site:
                    modifications["[NTerm_mods]"].append(mod_name)
                elif "C-term" in mod_site:
                    modifications["[CTerm_mods]"].append(mod_name)
                else:
                    modifications["[mods]"].append(mod_name)

    with open(param_path, "r", encoding="utf-8-sig") as param_file:
        lines = param_file.readlines()
    for section_name, prefix in (
        ("[mods]", "mod"),
        ("[NTerm_mods]", "Nmod"),
        ("[CTerm_mods]", "Cmod"),
    ):
        lines = _replace_mod_entries(
            lines, section_name, prefix, modifications[section_name]
        )
    with open(param_path, "w", encoding="utf-8", newline="") as param_file:
        param_file.writelines(lines)


def _read_modification_ini_names(modification_ini_path):
    known_names = set()
    with open(modification_ini_path, "r", encoding="utf-8-sig") as ini_file:
        for line in ini_file:
            line = line.strip()
            if (
                not line
                or line.startswith(("@", "#", "name"))
                or "=" not in line
            ):
                continue
            known_names.add(line.split("=", 1)[0].strip())
    return known_names


def fill_missing_mod_masses(param_path, modification_ini_path):
    known_names = _read_modification_ini_names(modification_ini_path)
    target_sections = {"[mods]", "[NTerm_mods]", "[CTerm_mods]"}
    current_section = None
    changed = False

    with open(param_path, "r", encoding="utf-8-sig") as param_file:
        lines = param_file.readlines()

    new_lines = []
    for line in lines:
        stripped = line.strip()
        if stripped.startswith("[") and stripped.endswith("]"):
            current_section = stripped
            new_lines.append(line)
            continue
        if (
            current_section not in target_sections
            or not stripped
            or stripped.startswith("#")
            or "=" not in stripped
        ):
            new_lines.append(line)
            continue

        key, value = stripped.split("=", 1)
        mod_name = value.split("\t", 1)[0].strip()
        if "\t" in value or mod_name in known_names:
            new_lines.append(line)
            continue

        while True:
            mass = input(f"Please enter the mass of modification {mod_name}: ").strip()
            try:
                float(mass)
                break
            except ValueError:
                print("Mass must be a number. Please try again.")
        line_ending = "\n" if line.endswith("\n") else ""
        new_lines.append(f"{key}={mod_name}\t{mass}{line_ending}")
        changed = True

    if changed:
        with open(param_path, "w", encoding="utf-8", newline="") as param_file:
            param_file.writelines(new_lines)


def read_summary(summarypath):
    m_str_lst = []
    m_mass_lst = []
    with open(summarypath, "r") as rsp:
        slines = rsp.readlines()
        for line in slines[1:]:
            line_lst = line.split("\t")
            m_str = line_lst[1]
            m_mass = line_lst[2].split(" ")[0]
            m_str_lst.append(m_str)
            m_mass_lst.append(m_mass)
    if len(m_str_lst) <= 3:
        return m_str_lst, m_mass_lst
    else:
        return m_str_lst[:3], m_mass_lst[:3]


def get_modstr(i_mod, is_mod, m_str_lst, m_mass_lst):
    if is_mod:
        mod = ""
        for i in range(len(m_str_lst)):
            for aa in aa_lst:
                m_str_ = m_str_lst[i] + "[" + aa + "]"
                mod += "mod" + str(i_mod) + "=" + m_str_ + "\t" + m_mass_lst[i] + "\n"
                i_mod += 1
        return mod
    else:
        modname = ""
        for i in range(len(m_str_lst)):
            for aa in aa_lst:
                m_str_ = m_str_lst[i] + "[" + aa + "]"
                modname += "ModName" + str(i_mod) + "=" + m_str_ + "\n"
                i_mod += 1
        return modname


def writepChemparam(parampath, wparampath, summarypath):
    m_str_lst, m_mass_lst = read_summary(summarypath)
    modsflag = False
    scoremodsflag = False
    modscount = 0
    with open(parampath, "r") as rpp:
        rlines = rpp.readlines()
        with open(wparampath, "w") as wpp:
            for line in rlines:
                if line.startswith("[mods]"):
                    modsflag = True
                if modsflag:
                    if line != "\n":
                        modscount += 1
                    else:
                        modstr = get_modstr(modscount, True, m_str_lst, m_mass_lst)
                        wpp.write(modstr)
                        modsflag = False
                        modscount = 0

                if line.startswith("[scoreModName]"):
                    scoremodsflag = True
                if scoremodsflag:
                    if line != "\n":
                        modscount += 1
                    else:
                        modstr = get_modstr(modscount, False, m_str_lst, m_mass_lst)
                        wpp.write(modstr)
                        scoremodsflag = False
                        modscount = 0
                wpp.write(line)


aa_lst = "ACDEFGHKLMNPQRSTVWY"


def read_embdict(summarypath):
    m_str_lst, m_mass_lst = read_summary(summarypath)
    modict = {}
    maxmodindex = 0
    with open(resource_path(r"pPred-MS2RT\emb_dict.txt"), "r") as emb:
        lines = emb.readlines()
        for line in lines:
            line = line.strip()
            modname = line.split("=")[0]
            modindex = line.split("=")[1]
            modict[modname] = modindex
            maxmodindex = int(modindex)
    for m_str in m_str_lst:
        if m_str not in modict.keys():
            with open(resource_path(r"pPred-MS2RT\emb_dict.txt"), "a") as wemb:
                for aa in aa_lst:
                    maxmodindex = maxmodindex + 1
                    xstr = m_str + "[" + aa + "]=" + str(maxmodindex)
                    wemb.write(xstr)
