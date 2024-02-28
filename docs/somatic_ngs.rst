.. _somatic_ngs:

-------------------
Mapping and friends
-------------------

To reliably identify somatic variants in low-purity, degraded samples, or when the focus is on sub-clonality, the quality of the mapping can become a limiting factor.
To make the best use of the data, it is sometimes necessary to consider adapter trimming & base quality recalibration.

Adapter trimming is currently implemented as a separate pipeline step, while the handling of UMIs (or Molecular BarCodes MBCs) and Base Quality Score Re-calibration (BQSR) have both been included in the ``ngs_mapping`` step.
Actually, UMIs & BQSR are included as sub-steps of the ``somatic`` mapping tool.
This "tool" is a placeholder for a collection of operations required to perform mapping & BQSR in presence of UMIs.

.. note:: 
    
    Adapter trimming should be included in this *meta-tool*, and eventually the whole operations carried out by the ``somatic`` tool should be done as the whole step.
    This is currently incompatible with the design of the ``ngs_mapping`` step, for efficiency reasons.

Adapter trimming
================

The ``ngs_mapping`` step provides a fast tool to trim adapters on-the-fly. 
However, there are cases where this fast but limited tool in not sufficient.
The ``adapter_trimming`` step provides 2 advanced tools (``bbduk`` & ``fastp``) for this task.

It is important to note that, unlike other pipeline steps, ``adapter_trimming`` produces ``fastq`` files.
This means that any subsequent pipeline step relying on ``fastq`` files for its input should be using the ``adapter_trimming`` output, not the file found in the standard way.

The configuration option ``path_link_in`` available for the ``ngs_mapping``, ``hla_typing`` ``ngs_data_qc`` & ``somatic_gene_fusion_calling`` steps must be used after ``adapter_trimming``.

The configuration snippet would then be similar to:

.. code-block:: yaml

    step_config:
      adapter_trimming:
        tools: [bbduk]
    
      ngs_mapping:
        path_link_in: <Absolute path to project folder>/adapter_trimming/output/bbduk
        ...

``path_link_in`` is substituted to the ``search_path`` entries from the ``data_sets`` sections. 
The pattern matching remains unchainged, as the ``adapter_trimming`` step does **not** rename any of the ``fastq`` files.

Barcodes & UMIs
===============

The handling of barcodes in generally done in 3 distinct opreations:

1. The barcodes (typically stored as the first few bases of the read) are clipped from the read sequencce and added to the read name or as comments to the description on the sequence identifer line of the read record.
2. Mapping is carried out as normal, but depending on the downstream tool, the MBC sequences & their qualities must be added as tags in the output ``bam`` file.
3. Aligned reads and the MBC sequences are used for de-duplication. This operation sometimes result in reads longer than actual reads.

The mapping tool can be selected from the mapping tools available for DNA.
However, the pipeline has **not** been tested with long reads mappers, only with mappers for short reads (Illumina), *i.e.* ``bwa`` & ``bwa-mem2``.
The arguments for the mapper are taken from the corresponding section in the configuration file.

.. note::

    Currently, only the `AGeNT <https://www.agilent.com/en/product/next-generation-sequencing/hybridization-based-next-generation-sequencing-ngs/ngs-software/agent-232879>`_ software from Agilent is implemented.
    Eventually, this commercial, non-free software should be replaced by open-source alternatives, for example `UMI-tools <https://umi-tools.readthedocs.io/en/latest/index.html>`_

An example of the configuration required for the ``somatic`` tool would be:

.. code-block:: yaml

    ngs_mapping:
      tools:
        dna: [mbcs]
      somatic:
        mapping_tool: bwa_mem2
        barcode_tool: agent         # Only agent is currently implemented
        use_barcodes: true
        recalibrate: true
      bwa_mem2:
        path_index: <path_to_bwa-mem2 indices>
        trim_adapter: false
        mark_duplicates: true
        split_as_secondary: true
        extra_args: ["-C"]          # Use ["-C"] when UMI/MBC are present, and processed with AGeNT, otherwise [""]
      agent:
        prepare:
          path: <path to AGeNT trimmer software>
          lib_prep_type: v2         # Check AGeNT documentation, must be one of "halo", "hs", "xt", "v2", "qxt"
          extra_args:
          - "-polyG 8"              # Check AGeNT documentation, trimming polyG tails
          - "-minFractionRead 50"   # Check AGeNT documentation, ignore heavily trimmed reads
        mark_duplicates:
          path: <path to AGeNT creak software>
          path_baits: <path to baits>
          consensus_mode: HYBRID    # Check AGeNT documentation, must be one of "SINGLE", "DUPLEX"m "HYBRID"
          input_filter_args:        # Check AGeNT documentation, input read filters on mapping & base qualities
          - "-mm 13"
          - "-mr 13"
          - "-mq 30"
          consensus_filter_args: [] # Check AGeNT documentation, filtering over consensus
          extra_args: []            # Check AGeNT documentation, extra arguments
      bqsr:
        common_variants: <path to common germline variants>  # For example small_exac_common_3.vcf from the GATK bucket

A summary of the ``AGeNT`` documentation is available on the cluster (``/fast/work/groups/cubi/projects/biotools/AGeNT/AGeNT ReadMe.pdf``).

.. note::

    The meta-tool is named ``somatic`` in the main description, but the mapper prefix in file name is ``mbcs``.
    Both choices are poor, and will be eventually changed.

Base Quality Score Re-calibration
=================================

`BQSR should generally be applied <https://gatk.broadinstitute.org/hc/en-us/articles/360035890531-Base-Quality-Score-Recalibration-BQSR>`_.
It is only implemented as part of the ``somatic/mbcs`` tool. 
However, it can be appiled even in the absence of UMIs or barcodes, just set the configuration option ``use_barcodes`` to ``false``.

