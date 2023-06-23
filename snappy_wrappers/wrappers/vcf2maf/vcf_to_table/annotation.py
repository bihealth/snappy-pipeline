import re
import typing

import vcfpy

import exceptions


class Annotation:
    """Process jannovar ANN or vep CSQ INFO entries

    The configuration file must contain the INFO identifier (typically ANN or CSQ), and
    regular expression to extract of column headers, and to split the annotation items.
    For example, the configuration for annotations created by jannovar would be:
    annotations:
        id: "ANN"
        extract: "^Functional annotations:'(.+)'$"
        split: "\\|"

    The vcf header must be first parsed, to extract the column headers of the various
    annotated fields. Then, the values of the records can be obtained.

    When there are multiple alternative alleles, annotations might be specific to a
    particular one. In that case, there must be one column specifying the alternative
    allele in the annotation. This information can be provided in the configuration by:
    annotations:
        id: "CSQ"
        allele: "Allele"
    """

    def __init__(
        self,
        annotation_id: str,
        extract: str,
        split: str,
        allele_column: str = None,
        minimize: bool = False,
    ):
        """Create annotation object

        :param annotation_id: annotation ID in the vcf INFO fields
        :param extract: regular expression string to extract column headers from the Description field
        :param split: regular expression string to split the extracted contents of the Description field into column headers
        :param allele_column: annotation column header to extract the alt allele (or None, when this information isn't present)
        :param minimize: minimization of variants to patch a VEP bug
        """
        if annotation_id is None or annotation_id == "":
            raise exceptions.MissingValue("Missing annotation ID")
        if not extract or not split:
            raise exceptions.MissingValue(
                'Missing patterns for extraction "{}" or for split "{}"'.format(extract, split)
            )
        self.annotation_id = annotation_id
        self.allele_column = allele_column
        self.columns = []
        self.extract = re.compile(extract)
        self.split = re.compile(split)
        self.minimize = minimize
        self.minimizer = None

    def parseHeader(self, header: vcfpy.Header):
        """Parse vcf header to extract annotation column headers

        :param header: vcf header as vcfpy.Header object
        """
        ann_header = header.get_info_field_info(self.annotation_id)
        if not ann_header:
            raise exceptions.MissingValue(
                'Annotation "{}" cannot be found in vcf INFO fields'.format(self.annotation_id)
            )
        m = self.extract.match(ann_header.description)
        if m:
            self.columns = self.split.split(m.group(1))
        else:
            raise exceptions.IllegalValue(
                'Column headers can\'t be extracted from description "{}" by pattern "{}"'.format(
                    ann_header.description, self.extract
                )
            )
        if self.allele_column:
            if self.allele_column not in self.columns:
                raise exceptions.MissingValue(
                    'Allele column "{}" not in annotation "{}" columns {}'.format(
                        self.allele_column, self.annotation_id, self.columns
                    )
                )

    def set_minimizer(self, record: vcfpy.Record):
        """Create a minimizer object for alternate alleles

        ENSEMBL VEP annotates vcf records with minimized alternative allele sequence,
        even without [`--minimal` option](https://www.ensembl.org/info/docs/tools/vep/script/vep_options.html).
        This undocumented feature is discussed in a [closed issue on github](https://github.com/Ensembl/ensembl-vep/issues/1264)

        To ensure that all variants are correctly processed, a minimizer must be created,
        to match the sequences of the input and minimised variants.
        This minimizer must be created for each record.

        :param record: vcf record
        """
        minimisation = record.REF
        for alt in record.ALT:
            a = alt.value
            i = 0
            while i < len(minimisation) and i < len(a) and minimisation[i] == a[i]:
                i += 1
            if i < len(minimisation):
                if i == 0:
                    minimisation = None
                    break
                else:
                    minimisation = minimisation[:i]
        if minimisation:
            self.minimizer = {}
            n = len(minimisation)
            for alt in record.ALT:
                if len(alt.value) == n:
                    self.minimizer["-"] = alt.value
                else:
                    self.minimizer[alt.value[n:]] = alt.value
        else:
            self.minimizer = None

    def getAnnotations(self, record: vcfpy.Record) -> typing.Dict[str, typing.List[typing.Any]]:
        """Extract annotations from a vcf record

        Depending on the choice of annotation software, there can be multiple annotations for a single variant.
        The annotations may also depend on the alternative allele (when more than one are present).
        The annotations are returned as a dict with the alternative alleles are keys.
        If the annotations are allele-dependent, they are attributed to the allele, otherwise
        all annotations are made available to all alleles.
        For each allele, the annotations are given as a dict, with the column name as key,
        and the list of annotations for that column as value.

        :param record: vcf record as vcfpy.Record object
        :return: dict of of dict of annotation lists
        """
        if record is None:
            raise exceptions.MissingValue("Missing vcf record")
        ann = record.INFO[self.annotation_id]
        if not isinstance(ann, list):
            ann = [ann]
        alts = [alt.value for alt in record.ALT]
        values = {}
        for theColumn in self.columns:
            values[theColumn] = {alt: None for alt in alts}
        for a in ann:
            one_ann = self._get_one_annotation(a)
            storeAlts = self.columns
            if self.allele_column and one_ann[self.allele_column] in alts:
                storeAlts = [one_ann[self.allele_column]]
            for k, v in one_ann.items():
                for alt in storeAlts:
                    if values[k][alt] is None:
                        values[k][alt] = []
                    values[k][alt].append(v)

        return values

    def _get_one_annotation(self, annotation: str) -> typing.Dict[str, str]:
        values = self.split.split(annotation)
        assert len(values) == len(self.columns)
        values = [x if x else None for x in values]
        result = dict(zip(self.columns, values))
        if self.allele_column and self.minimizer and result[self.allele_column] in self.minimizer:
            result[self.allele_column] = self.minimizer[result[self.allele_column]]
        return result
