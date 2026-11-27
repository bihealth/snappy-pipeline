{% if snakemake.config.get("__rulegraph__") %}
The rules of this run as a metro map: each line follows the outputs of one rule to the rules that read them.
The jobs that link the results and the final target are left out.

.. embeddedimage:: {{ snakemake.config["__rulegraph__"] }}
{% else %}
A full ``snappy run`` draws the rules of the run as a metro map, which the report shows here.
{% endif %}
