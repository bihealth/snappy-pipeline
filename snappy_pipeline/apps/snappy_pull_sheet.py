"""Pull the sample sheet of a SODAR project through its API, as a plain TSV sheet."""

import sys

import requests

#: The URL template to use; it holds the API key, so it is never printed.
URL_TPL = "https://%(sodar_host)s/samplesheets/api/remote/get/%(project_uuid)s/%(api_key)s"

#: Seconds to wait for the SODAR server to connect and to send data.
TIMEOUT = 60

#: The columns of the written sheet
HEADER = (
    "family_id\tdonor_name\tfather\tmother\tsex\tis_affected\tsample_name\tlibrary_name"
    "\tlibrary_type\tfolder_name\tbatch_no\tseq_platform\tlibrary_kit\tproject_name\tproject_uuid"
)

#: Mapping from ISA-tab sex to the sheet's sex.
MAPPING_SEX = {"female": "female", "male": "male", "unknown": "unknown", None: ""}

#: Mapping from disease status to the sheet's affected status.
MAPPING_STATUS = {
    "affected": "true",
    "carrier": "true",
    "unaffected": "false",
    "unknown": "",
    None: "",
}


def strip(x):
    if hasattr(x, "strip"):
        return x.strip()
    else:
        return x


def run(args):  # noqa: C901
    # Query investigation JSON from API.
    print(
        "Fetching the sample sheet of project %s from %s" % (args.project_uuid, args.sodar_host),
        file=sys.stderr,
    )
    try:
        r = requests.get(URL_TPL % vars(args), timeout=TIMEOUT)
        r.raise_for_status()
    except requests.RequestException as e:
        # The messages of requests name the URL, and with it the API key
        status = f" {e.response.status_code}" if e.response is not None else ""
        raise SystemExit(f"Fetching the sample sheet failed: {type(e).__name__}{status}") from None
    all_data = r.json()
    if len(all_data["studies"]) > 1:
        raise Exception("More than one study found!")

    # Parse out study data.
    study = list(all_data["studies"].values())[0]
    study_infos = study["study"]
    study_top = study_infos["top_header"]
    n_source = study_top[0]["colspan"]
    n_extraction = study_top[1]["colspan"]
    cols_source = study_infos["field_header"][:n_source]
    cols_extraction = study_infos["field_header"][n_source : n_source + n_extraction]
    cols_sample = study_infos["field_header"][n_source + n_extraction :]
    names_source = [x["value"] for x in cols_source]
    names_extraction = [x["value"] for x in cols_extraction]
    names_sample = [x["value"] for x in cols_sample]
    table = study_infos["table_data"]

    # Build study info map.
    study_map = {}
    for row in table:
        # Assign fields to table.
        dict_source = dict(zip(names_source, [strip(x["value"]) for x in row[:n_source]]))
        dict_extraction = dict(
            zip(names_extraction, [x["value"] for x in row[n_source : n_source + n_extraction]])
        )
        dict_sample = dict(
            zip(names_sample, [strip(x["value"]) for x in row[n_source + n_extraction :]])
        )
        # Extend study_map.
        study_map[dict_source["Name"]] = {
            "Source": dict_source,
            "Extraction": dict_extraction,
            "Sample": dict_sample,
        }

    # Parse out the assay data.
    #
    # NB: We're not completely cleanly decomposing the information and, e.g., overwrite
    # the "Extract name" keys here...
    if len(study["assays"]) > 1:
        raise Exception("More than one assay found!")
    assay = list(study["assays"].values())[0]
    top_columns = [(x["value"], x["colspan"]) for x in assay["top_header"]]
    columns = []
    offset = 0
    for type_, colspan in top_columns:
        columns.append(
            {
                "type": type_,
                "columns": [x["value"] for x in assay["field_header"][offset : offset + colspan]],
            }
        )
        offset += colspan
    assay_map = {}
    for row in assay["table_data"]:
        offset = 0
        name = row[0]["value"].strip()
        for column in columns:
            colspan = len(column["columns"])
            values = {
                "type": column["type"],
                **dict(
                    zip(column["columns"], [x["value"] for x in row[offset : offset + colspan]])
                ),
            }
            type_ = column["type"]
            if type_ == "Process":
                type_ = values["Protocol"]
            assay_map.setdefault(name, {})[type_] = values
            offset += colspan

    # Generate the resulting sample sheet.
    print(HEADER, file=args.output)
    for source, info in study_map.items():
        if source not in assay_map:
            print("INFO: source %s does not have an assay." % source, file=sys.stderr)
            dict_lib = {}  # a donor without a library
            proc_lib = {}
        else:
            for outer_key in ("Extract Name", "Library Name"):
                if outer_key in assay_map[source]:
                    dict_lib = assay_map[source][outer_key]
                    for key in assay_map[source]:
                        if key.startswith("Library construction"):
                            proc_lib = assay_map[source][key]
                            break
                    else:
                        proc_lib = {}
        dict_source = info["Source"]
        # FIXME: remove hack
        # HACK: ignore if looks like artifact
        if "Folder Name" in dict_lib:
            library_name = dict_lib["Name"]
            library_type = library_name.split("-")[-1][:-1]  # hack to get library type
            folder = dict_lib["Folder Name"]
        else:
            library_name = library_type = folder = ""
        # TODO: find better way of accessing sequencing process
        seq_platform = "Illumina" if library_name else ""
        for the_dict in assay_map.get(source, {}).values():
            if the_dict.get("Platform") == "PACBIO_SMRT":
                seq_platform = "PacBio"
        library_kit = proc_lib.get("Library Kit") or ""
        if args.library_types and library_type and library_type not in args.library_types:
            print(
                "Skipping %s not in library types %s" % (dict_source["Name"], args.library_types),
                file=sys.stderr,
            )
            continue
        # ENDOF HACK
        parents = [dict_source.get(key) or "" for key in ("Father", "Mother")]
        row = [
            dict_source["Family"],
            dict_source["Name"],
            *("" if parent == "0" else parent for parent in parents),
            MAPPING_SEX.get((dict_source.get("Sex") or "").lower(), "unknown"),
            MAPPING_STATUS.get((dict_source.get("Disease Status") or "").lower(), ""),
            info["Sample"].get("Name", "") if library_name else "",
            library_name,
            library_type,
            folder,
            dict_source.get("Batch") or "",  # the batch of SODAR verbatim: an id, not a number
            seq_platform,
            library_kit,
            args.project_name,
            args.project_uuid,
        ]
        print("\t".join(row), file=args.output)
