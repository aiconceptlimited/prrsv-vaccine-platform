# VAXINTAIC Scientific Status

## Scope

VAXINTAIC is a computational framework for multi-stage PRRSV vaccine
candidate design and in-silico prioritization.

The framework integrates sequence analysis, epitope screening, candidate
ranking, mRNA construct design, host-specific codon adaptation,
computational immunogenicity scoring, diversity analysis, and modelled
delivery/prioritization metrics.

## Evidence status

Repository outputs are computational results unless explicitly identified
as publication-derived materials.

Computational scores must not be interpreted as experimental measurements
of vaccine efficacy, protection, immunogenicity, safety, or nanoparticle
delivery.

## Nanoparticle modelling

The nanoparticle/LNP component is a computational modelling component.
Its training procedure uses simulated training data.

Its outputs therefore do not constitute experimental evidence of lipid
nanoparticle encapsulation, delivery, vaccine efficacy, or protection.

## Computational prioritization

The candidate model produces in-silico prioritization metrics.

These metrics are not calibrated probabilities of vaccine efficacy and
should not be presented as experimentally validated predictions of
protective efficacy.

## Immunogenicity

Immunogenicity scoring is computational and sequence-derived. It should
not be interpreted as experimentally measured immunogenicity.

## Diversity

The current diversity metric is a dominance-complement metric:

    diversity = 1 - dominant_frequency

It is not Shannon entropy, Simpson diversity, or another standardized
diversity index.

## Codon adaptation

Codon adaptation uses a documented Sus scrofa host-specific codon-use
reference. CAI is only reported when the required reference is available.

## Publication materials

The VAXINTAIC TVJ publication figures and associated submission materials
are retained as publication-derived materials.

Their accompanying validation records document the source records,
sequence checks, construct checks, and figure-generation status.

## Reproducibility

The authoritative pipeline runner records run-level provenance,
including execution metadata and artifact information.

Runtime credentials, private databases, local virtual environments,
machine-specific paths, logs, backups, and runtime history are not part
of the public source release.

## Interpretation

VAXINTAIC supports computational vaccine candidate design and
prioritization.

Experimental studies remain necessary to establish biological activity,
immunogenicity, delivery performance, safety, and protective efficacy.
