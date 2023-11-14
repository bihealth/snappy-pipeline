import logging
import re
import typing

from annotation import Annotation
import exceptions
import vcfpy


class VcfParser:
    pattern = re.compile("^([^\[\]\s]+)(\[([0-9]+|REF)\])?$")  # noqa: W605

    def __init__(
        self,
        config: typing.Dict[str, typing.Any],
        header: vcfpy.Header,
        samples: typing.Dict[str, str],
    ):
        """Create the parser for column requests

        The parser is configured using the configuration yaml file, in particular the "annotation" dict,
        and the info & format list of requested IDs.
        The latter is used to give the user control of the set of INFO & FORMAT IDs that can be used in parsing.
        When these lists are missing from the configuration, all the INFO & FORMAT IDs present in the VCF can be used.

        The annotation defines the extraction pattern for parsing the annotation field description in the VCF header.
        The description field generally consists of a prefix, followed by the column names of the annotations.
        The `extract` regular expression must define the column names as a single string in the first group,
        and the `split` regular expression is used to separate the annotation field (in the header or in the records)
        into individual values.

        The sample names must be provided as a dict with `tumor` and `normal` keys (the latter being optional).
        This provides a mapping between logical labels (`tumor` & `normal`) and the sample column names.

        :param config: configuration dict
        :param header: VCF header as vcfpy.Header object
        :param samples: samples mapping dict
        """
        if samples["normal"]:
            if not samples["normal"] in header.samples.names:
                raise exceptions.IllegalValue(
                    'Normal sample "{}" not found in vcf header'.format(samples["normal"])
                )
        if samples["tumor"]:
            if not samples["tumor"] in header.samples.names:
                raise exceptions.IllegalValue(
                    'Tumor sample "{}" not found in vcf header'.format(samples["tumor"])
                )
        self.samples = samples

        self.annotation = None
        self.infos = VcfParser._get_id_list(config.get("info", None), header, "INFO")
        self.formats = VcfParser._get_id_list(config.get("format", None), header, "FORMAT")

        if (
            "annotation" in config
            and isinstance(config["annotation"], dict)
            and set(["id", "extract", "split"]) <= set(config["annotation"].keys())
        ):
            self.annotation = Annotation(
                annotation_id=config["annotation"]["id"],
                allele_column=config["annotation"]["allele"]
                if "allele" in config["annotation"]
                else None,
                extract=config["annotation"]["extract"],
                split=config["annotation"]["split"],
            )
            self.annotation.parseHeader(header)

        self.constants = {}

    def addConstant(self, key, value):
        """Add an key, value pair to the parser set of constant values

        :param key: name of the constant column
        :param value: value of all elements in the constant column
        """
        if key in self.constants:
            logging.warning(
                'Value of constant "{}" is replaced from {} to {}'.format(
                    key, self.constants[key], value
                )
            )
        self.constants[key] = value

    @staticmethod
    def _get_id_list(
        request, header, line_type=""
    ) -> typing.Dict[str, vcfpy.header.CompoundHeaderLine]:
        """Creates a dict of FORMAT or INFO lines

        The user can subset the ID found in the vcf using the configuration file.
        If a list of requested IDs is provided, these IDs will be checked for presence in the VCF.
        If any requested ID are missing from the VCF, an exception is raised.
        If there is no user request, all IDs present in the VCF header are returned.

        :param request: list of ids, possibly empty or None
        :param header: VCF header as vcfpy object
        :param line_type: either INFO or FORMAT. Other values trigger an exception
        :return: dict with ID as keys, and header lines as values (vcfpy.CompoundHeaderLine)
        """
        if line_type == "INFO":
            ids = header.info_ids()
            value_fct = header.get_info_field_info
        elif line_type == "FORMAT":
            ids = header.format_ids()
            value_fct = header.get_format_field_info
        else:
            raise exceptions.IllegalValue('Illegal VCF header line "{}"'.format(line_type))
        id_list = {x: value_fct(x) for x in ids}
        if request:
            request = set(request)
            in_header = set(id_list.keys())
            missing = request - in_header
            if len(missing) > 0:
                raise exceptions.MissingValue(
                    'Requested {} id(s) not found in vcf header: "{}"'.format(
                        '", "'.join(missing), '", "'.join(in_header)
                    )
                )
            unused = in_header - request
            for field in unused:
                del id_list[field]
        return id_list

    def subset(self, record: vcfpy.Record, index=None) -> vcfpy.Record:
        """Subset a VCF record to retain on one alt allele

        A new VCF record is created, to represent the ith alternative allele of the variant.
        The allele selection is provided by the index (0 is the first alt allele).
        When the index is None, the original VCF record is returned, with all the alternative alleles.
        When the index is greater than the number of alt alleles found for the variant,
        an empty record is returned (None).

        :param record: the VCF record with multiple variants
        :param index: the index of the selected alt allele (starting from 0)
        :return: the variant with only the selected alt allele, or an empty record when the ith allele is missing, or the original record if the index is missing
        """
        if index is None:
            return record
        if index >= len(record.ALT):
            return None
        infos = {}
        for ID in record.INFO.keys():
            assert ID in self.infos
            header = self.infos[ID]
            if header.number == vcfpy.header.HEADER_NUMBER_REF:
                infos[ID] = [record.INFO[ID][0], record.INFO[ID][1 + index]]
            elif header.number == vcfpy.header.HEADER_NUMBER_ALLELES:
                infos[ID] = [record.INFO[ID][index]]
            else:
                infos[ID] = record.INFO[ID]

        formats = []
        for call in record.calls:
            data = {}
            for ID in record.FORMAT:
                assert ID in self.formats
                fmt = self.formats[ID]
                if fmt.number == vcfpy.header.HEADER_NUMBER_REF:
                    data[ID] = [call.data[ID][0], call.data[ID][1 + index]]
                elif fmt.number == vcfpy.header.HEADER_NUMBER_ALLELES:
                    data[ID] = [call.data[ID][index]]
                else:
                    data[ID] = call.data[ID]
            formats.append(vcfpy.record.Call(call.sample, data))

        return vcfpy.record.Record(
            record.CHROM,
            record.POS,
            record.ID,
            record.REF,
            [record.ALT[index]],
            record.QUAL,
            record.FILTER,
            infos,
            record.FORMAT,
            formats,
        )

    def getInputValue(
        self,
        record: vcfpy.Record,
        theSet: str,
        theColumn: str,
        theSample: str = None,
        default: str = None,
    ):
        """Extract the requested value(s) from a VCF record

        The request format depends on the set of values the data must be extracted from.
        The available sets are:

        - `constant`: constant values stored in the parser (typically the genome release, or the organisation name).
           The request must be the constant name.
        - `fixed`: the fixed columns of the VCF. The request must be one of `CHROM`, `POS`, ÌD`, `QUAL` or `FILTER`.
        - `variant`: the reference or alternative variant. The request must be `REF` or `ALT`, or `ALT[<index>]`.
          If the request is `ALT`, then all alt variants will be returned.
          If the request is `ALT[<index>]`, then only the ith variant is returned, and if there is no ith variant, None is returned.
        - `info`: an INFO field, by ID. The request must be `<the ID>`, `<the ID>[<index>]` or `<the ID>[REF]`.
          The meaning of the index (either a positive integer, starting from 0, or `REF`) depends on the
          INFO field numbering scheme, as defined in the VCF header.
          When the number of values depends on the number of alt alleles (either `R` or `A` in VCF Number parlance),
          only one value can be returned. If the value list contains the reference (`R` number of values),
          the value for the reference can be obtained using the `REF` index.
          Alternatively, when the number of values does not depend on the number of alt alleles,
          then the indexing can be used. When it isn't used, all values are returned.
        - `format`: a `FORMAT` field, by ID. The request format is identical to the format used for `INFO` fields,
          but the logical sample id (either `tumor` or `normal`) must be provided.
        - `annotation`: an annotation field, by column name. The request format is identical to the format used for `INFO` fields.
          When annotations have been configured to extract the allele from the annotation columns, then
          the annotations are assumed to have a `A` number, so they are attached to alt alleles. No indexing is possible.
          If it is not the case, the annotations can have any number of values. Indexing by number is possible.

        The returned values are either a list (of length 1 for `constant` and `fixed` sets), or a dict of lists
        when the requested values depend on the alt allele. In that case, the lists lengths are always 1 for
        `INFO` and `FORMAT` fields, but can be longer for annotations (as there can be multiple annotations per alt allele).

        :param record: the VCF record the data is extracted from, as a vcfpy.Record object
        :param theSet: one of `constant`, `fixed`, `variant`, `info`, `format` or `annotation`, to describe the type of data for extraction
        :param theColumn: the column ID, possibly with indexing
        :param theSample: sample logical label (`tumor` or `normal`), for `format` data extraction
        :param action: action to be taken for missing values. One of `ignore`, `skip`, `stop`, `default`.
        :return: the parsed values, as a dict for alt-allele dependent values, or a list otherwise, or None when absent
        """
        if record is None:
            raise exceptions.MissingValue("No record provided")
        if not theSet or not theColumn:
            raise exceptions.MissingValue(
                'Missing set or column information for record "{}"'.format(record)
            )
        if theSet == "format":
            if not theSample:
                raise exceptions.MissingValue(
                    'Missing sample information for column "{}" and record "{}"'.format(
                        theColumn, record
                    )
                )
            if theSample not in self.samples:
                raise exceptions.IllegalValue('Unknown sample "{}"'.format(theSample))
            theSample = self.samples[theSample]

        value = None
        if theSet == "constant":
            if theColumn in self.constants:
                value = self.constants.get(theColumn, None)
                if value is not None:
                    value = VcfParser._expandUnique(value, record.ALT)
        elif theSet == "fixed":
            value = VcfParser._expandUnique(VcfParser._getFixed(record, theColumn), record.ALT)
        elif theSet == "variant":
            value = VcfParser._getVariant(record, theColumn)
        elif theSet == "info":
            numbers = {x.id: x.number for x in self.infos.values()}
            infos = VcfParser._expandDict(record.INFO, record.ALT, numbers)
            value = self._getRequest(theColumn, infos, numbers)
        elif theSet == "format":
            if theSample:
                numbers = {x.id: x.number for x in self.formats.values()}
                formats = VcfParser._expandDict(
                    record.call_for_sample[theSample].data, record.ALT, numbers
                )
                value = self._getRequest(theColumn, formats, numbers)
        elif theSet == "annotation":
            numbers = {x: vcfpy.header.HEADER_NUMBER_UNBOUNDED for x in self.annotation.columns}
            annotations = self.annotation.getAnnotations(record)
            value = self._getRequest(theColumn, annotations, numbers)
        else:
            raise exceptions.IllegalValue('Unknown set "{}"'.format(theSet))

        if default is not None:
            value = VcfParser._insert_default(value, default, [alt.value for alt in record.ALT])
        else:
            assert value is None or isinstance(value, dict)
            if isinstance(value, dict):
                for k, v in value.items():
                    assert v is None or isinstance(v, list)

        return value

    @staticmethod
    def _insert_default(value: typing.Dict[str, typing.Any], default: str, alts: typing.List[str]):
        if value is None:
            value = {alt: [default] for alt in alts}
        else:
            assert isinstance(value, dict)
            for alt in alts:
                if alt in value:
                    if value[alt] is None:
                        value[alt] = [default]
                else:
                    assert isinstance(value[alt], list)
                    for i in range(len(value[alt])):
                        if value[alt][i] is None:
                            value[alt][i] = default
        return value

    @staticmethod
    def _expandUnique(
        value, alts: typing.List[vcfpy.record.AltRecord]
    ) -> typing.Dict[str, typing.Any]:
        if value is None:
            return {alt.value: value for alt in alts}
        else:
            return {alt.value: [value] for alt in alts}

    @staticmethod
    def _expandDictElement(
        value, alts: typing.List[vcfpy.record.AltRecord], number
    ) -> typing.Dict[str, typing.Any]:
        if value is None:
            return None
        if not isinstance(value, list):
            return {alt.value: [value] for alt in alts}
        if number == vcfpy.header.HEADER_NUMBER_UNBOUNDED:
            return {alt.value: value for alt in alts}
        if number == vcfpy.header.HEADER_NUMBER_ALLELES:
            return {alts[i].value: [value[i]] for i in range(len(alts))}
        if number == vcfpy.header.HEADER_NUMBER_REF:
            result = {alts[i].value: [value[i + 1]] for i in range(len(alts))}
            result["__REF__"] = [value[0]]
            return result
        return {alt.value: value for alt in alts}

    @staticmethod
    def _expandDict(
        value, alts: typing.List[vcfpy.record.AltRecord], numbers
    ) -> typing.Dict[str, typing.Dict[str, typing.Any]]:
        if value is None:
            return None
        return {k: VcfParser._expandDictElement(v, alts, numbers[k]) for k, v in value.items()}

    @staticmethod
    def _getFixed(record: vcfpy.Record, theColumn: str) -> typing.List[typing.Any]:
        """Extract the fixed part of a vcf record

        The variant position (CHROM & POS), its identifier (ID), filter (FILTER) & quality (QUAL) values can be returned.
        The allowed column values are: CHROM, POS, ID, FILTER & QUAL.
        Note that multiple IDs & multiple filters are collapsed into a single string with the semicolumn (;) as separator

        :param record: vcf record as vcfpy.Record object
        :param theColumn: which fixed value must be returned
        :return: an array of length 1 with the requested fixed value
        """
        if theColumn == "CHROM":
            return record.CHROM
        elif theColumn == "POS":
            return record.POS
        elif theColumn == "ID":
            return ";".join(record.ID)
        elif theColumn == "QUAL":
            return record.QUAL
        elif theColumn == "FILTER":
            return ";".join(record.FILTER)
        else:
            raise exceptions.IllegalValue('Unknown fixed column "{}"'.format(theColumn))

    @staticmethod
    def _getVariant(record: vcfpy.Record, theColumn: str) -> typing.Dict[str, typing.List[str]]:
        """Extract the reference and alternative allele(s) from a vcf record

        A single reference & possibly multiple alternative alleles are returned as a list.
        The allowed column names are 'REF' for the reference, and
        'ALT' for the first alt allele,'ALT[<number>]' for the nth, or 'ALT[*]' for all the alt alleles.
        Note that 'ALT[<number>]' will return an empty list in the requested number is larger than the number of alternative alleles for the variant.

        :param record: vcf record as vcfpy.Record object
        :param theColumn: the type of allele list to be returned
        :return: a dict of lists of length 1 for the reference allele, or a possibly longer list for alternative alleles
        """
        if theColumn == "REF":
            return {alt.value: [record.REF] for alt in record.ALT}
        m = VcfParser.pattern.match(theColumn)
        if not m:
            raise exceptions.IllegalValue('Unknown variant column "{}"'.format(theColumn))
        groups = m.groups()
        if groups[0] != "ALT":
            raise exceptions.IllegalValue('Unknown variant column "{}"'.format(theColumn))
        if groups[2]:
            value = {alt.value: None for alt in record.ALT}
            n = int(groups[2])
            if n < len(record.ALT):
                value[record.ALT[n].value] = [record.ALT[n].value]
            return value
        else:
            return {alt.value: [alt.value] for alt in record.ALT}

    def _getRequest(
        self,
        request: str,
        values: typing.Any,  # A dict with alt alleles as keys and with data relevant to the allele as value
        numbers: typing.Dict[str, str],
    ) -> typing.Dict[str, typing.Any]:
        if values is None:
            return None

        m = VcfParser.pattern.match(request)
        if not m:
            raise exceptions.IllegalValue('Syntax error in column request "{}"'.format(request))

        groups = m.groups()
        theColumn = groups[0]
        if theColumn not in values.keys():
            return None
        result = values[theColumn]

        if groups[2]:
            if groups[2] == "REF":
                if numbers[theColumn] != vcfpy.header.HEADER_NUMBER_REF:
                    raise exceptions.IllegalValue(
                        'Requested data for reference allele in column "{}" which doesn\'t contain it'.format(
                            theColumn
                        )
                    )
                assert "__REF__" in result
                ref_data = result.pop("__REF__")
                for alt in result.keys():
                    result[alt] = ref_data
            else:
                n = int(groups[2])
                for alt, v in result.items():
                    if n < len(v) and v[n] is not None:
                        result[alt] = [v[n]]
                    else:
                        result[alt] = None

        return result
