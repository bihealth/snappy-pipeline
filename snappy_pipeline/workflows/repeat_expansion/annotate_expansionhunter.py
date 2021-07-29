# -*- coding: utf-8 -*
import json
import os

#: Dictionary with possible output sentences
SENTENCES = {
    "normal": "Normal - repeat {} observed: {}. Normal range: {}. Expansion range: {}.",
    "expanded": "Expanded - repeat {} observed: {}. Normal range: {}. Expansion range: {}.",
    "unclear": (
        "Unclear - repeat {} observed outside defined ranges ({}). "
        "Normal range: {}. Expansion range: {}."
    ),
    "no_annotation": "Undefined - no annotation found for repeat {}. Observed genotype: {}.",
    "no_data": "Undefined - genotype not defined in ExpansionHunter for repeat {}.",
}


class AnnotateExpansionHunter:
    """Class contain methods required to annotate ExpansionHunter results."""

    def __init__(self, eh_json, annotation_json, output_path):
        """Constructor.

        :param eh_json: Path to ExpansionHunter output file - JSON format.
        :type eh_json: str

        :param annotation_json: Path to repeat expansion annotation file - JSON format.
        :type annotation_json: str

        :param output_path: Path to output file.
        :type output_path: str
        """
        # Validate input
        if not (os.path.isfile(eh_json) and os.access(eh_json, os.R_OK)):
            raise ValueError("Cannot access ExpansionHunter results file: {0}".format(eh_json))
        if not (os.path.isfile(annotation_json) and os.access(annotation_json, os.R_OK)):
            msg = "Cannot access repeat expansion annotation file: {0}".format(annotation_json)
            raise ValueError(msg)
        directory_path = os.path.dirname(os.path.abspath(output_path))
        if not (os.access(directory_path, os.W_OK)):
            raise ValueError("No permission to write to directory {0}".format(directory_path))
        # Set variables
        self.eh_json = eh_json
        self.annotation_json = annotation_json
        self.output_path = output_path

    def run(self):
        """Performs class implementation."""
        # Initialise variable
        output = {}

        # Load EH result
        eh_dict = self.load_json(path=self.eh_json)
        simplified_out_dict = eh_dict.get("LocusResults")
        sample_parameters_dict = eh_dict.get("SampleParameters")
        # Load annotation file
        annotation_dict = self.load_json(path=self.annotation_json)

        # Annotate repeats
        annotated_repeats = self.annotate(
            output_dict=simplified_out_dict, annotation_dict=annotation_dict
        )

        # Custom out
        output["annotation"] = annotation_dict
        output["results"] = annotated_repeats
        output["raw"] = simplified_out_dict
        output["sample_parameters"] = sample_parameters_dict

        # Write to file
        self.write_json(path=self.output_path, data=output)

    def annotate(self, output_dict, annotation_dict):
        """Annotate ExpansionHunter output.

        :param output_dict: Dictionary with ExpansionHunter JSON output.
        :type output_dict: dict

        :param annotation_dict: Dictionary with Repeat Expansions annotations.
        :type annotation_dict: dict

        :return: Returns dictionary with explanation/annotation for each repeat.
        """
        # Initialise variable
        explanation_dict = {}

        # Iterate over results
        for repeat in output_dict:
            # Initialise ranges
            normal_start = None
            normal_end = None
            expansion_start = None
            expansion_end = None
            annotation = annotation_dict.get(repeat, None)
            if annotation:
                normal_start = annotation.get("normal_range_start")
                normal_end = annotation.get("normal_range_end")
                expansion_start = annotation.get("expansion_range_start")
                expansion_end = annotation.get("expansion_range_end")

            # Iterate over variants
            # Example: repeat 'ATXN7' has two variants, 'ATXN7' and 'ATXN7_GCC'
            for repeat_variant in output_dict[repeat]["Variants"]:
                genotype = output_dict[repeat]["Variants"][repeat_variant].get("Genotype", None)
                # Update repeat name to accommodate for variants
                if repeat_variant != repeat:
                    repeat_variant = "{rep} ({var})".format(rep=repeat, var=repeat_variant)
                sentence = self.explain_result(
                    repeat=repeat_variant,
                    genotype=genotype,
                    normal_start=normal_start,
                    normal_end=normal_end,
                    expansion_start=expansion_start,
                    expansion_end=expansion_end,
                )
                # Update out dict
                explanation_dict[repeat] = sentence

        # Return explanation dict
        return explanation_dict

    def explain_result(
        self,
        repeat=None,
        genotype=None,
        normal_start=None,
        normal_end=None,
        expansion_start=None,
        expansion_end=None,
    ):
        """Explain result.

        :param repeat: Repeat expansion name, commonly gene name (e.g., 'HTT' and 'FMR1').
        :type repeat: str

        :param genotype: Number of copies as outputted by ExpansionHunter.
        :type genotype: str

        :param normal_start: Start of normal range.
        :type normal_start: int

        :param normal_end: End of normal range.
        :type normal_end: int

        :param expansion_start: Start of expanded range.
        :type expansion_start: int

        :param expansion_end: End of expanded range.
        :type expansion_end: int

        :return: Return sentence explain the result (i.e., if normal, expanded or unclear) together
        with the ranges used.
        """
        # Define ranges
        range_template = "{s} - {e}"
        normal_range_str = range_template.format(s=normal_start, e=normal_end)
        if expansion_end:
            expansion_end_str = range_template.format(s=expansion_start, e=expansion_end)
        else:
            expansion_end_str = "{s}+".format(s=expansion_start)

        # Get max genotype
        value = self._get_max_genotype(genotype)

        # No annotation available
        if not all([normal_start, normal_end, expansion_start]):
            sentence = SENTENCES.get("no_annotation")
            sentence = sentence.format(repeat, genotype)

        # No data reported
        elif not genotype:
            sentence = SENTENCES.get("no_data")
            sentence = sentence.format(repeat)

        # Normal range
        elif normal_start <= value <= normal_end:
            sentence = SENTENCES.get("normal")
            sentence = sentence.format(repeat, genotype, normal_range_str, expansion_end_str)

        # Expanded range
        elif expansion_start <= value <= expansion_end:
            sentence = SENTENCES.get("expanded")
            sentence = sentence.format(repeat, genotype, normal_range_str, expansion_end_str)

        # Unclear
        else:
            sentence = SENTENCES.get("unclear")
            sentence = sentence.format(repeat, genotype, normal_range_str, expansion_end_str)

        # Return sentence
        return sentence

    @staticmethod
    def _get_max_genotype(genotype):
        """Get max genotype number.

        :param genotype: Genotype value as outputted by ExpansionHunter, e.g.: '24/25'.

        :return: Returns max count found in genotype.
        """
        try:
            return int(genotype)
        except ValueError:
            return max([int(value) for value in genotype.split("/")])
        except TypeError:
            return None

    @staticmethod
    def load_json(path):
        """Load JSON file.

        :param path: Path to JSON file.
        :type path: str

        :return: Dictionary with JSON file representation.
        """
        with open(path) as json_file:
            return json.load(json_file)

    @staticmethod
    def write_json(path, data):
        """Write JSON file.

        :param path: Path to JSON file.
        :type path: str

        :param data: JSON file data/content.
        :type data: dict
        """
        with open(path, "w") as json_file:
            json.dump(data, json_file)
