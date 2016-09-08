#! /usr/bin/env python

import os
import sys
import collections
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, SCRIPT_DIR)
import s00
import csv

import random
import argparse
from dendropy.model import reconcile
from dendropy.interop import seqgen
import dendropy

BPP_TEMPLATE = """\

          seed =  -1

       seqfile = {chars_filepath}
      Imapfile = {imap_filepath}
       outfile = {out_filepath}
      mcmcfile = {mcmc_filepath}

* speciesdelimitation = 0 * fixed species tree
 speciesdelimitation = 1 0 2.5    * speciesdelimitation algorithm0 and finetune(e)
*  speciesdelimitation = 1 1 2 1  * speciesdelimitation algorithm1 finetune (a m)
*           speciestree = 1

     speciesmodelprior = 1         * 0: uniform labeled histories; 1:uniform rooted trees

  species&tree = {num_species}  {species_labels}
                     {num_individuals_per_species}
                 {species_tree}

       usedata = 1    * 0: no data (prior); 1:seq like
         nloci = {num_loci}    * number of data sets in seqfile

     cleandata = 0    * remove sites with ambiguity data (1:yes, 0:no)?

    thetaprior = {theta_prior_a} {theta_prior_b}   # gamma(a, b) for theta
      tauprior = {tau_prior_a}   {tau_prior_b}       # gamma(a, b) for root tau & Dirichlet(a) for other tau's

*       finetune = 0: 5 0.0005 0.002  0.0005 0.5 0.2 1.0  # auto (0 or 1): finetune for GBtj, GBspr, theta, tau, mix, locusrate, seqerr

       finetune = 1: .01 .01 .01 .01 .01 .01 .01 .01  # auto (0 or 1): finetune for GBtj, GBspr, theta, tau, mix, locusrate, seqerr

         print = 1 0 0 0   * MCMC samples, locusrate, heredityscalars Genetrees
        burnin = 4000
      sampfreq = 2
       nsample = 10000
"""

def try_to_coerce_to_float(v):
    try:
        return float(v)
    except ValueError:
        if v == "None":
            return "NA"
        else:
            return v

def multinomial_sample(n, pvals):
    import numpy
    k = numpy.random.multinomial(n, pvals)
    return list(k)

def process_containing_tree_for_gene_samples(
        containing_tree,
        total_number_of_individuals=200,
        rng=None):
    n = total_number_of_individuals
    # pvals = [1.0/len(containing_tree.taxon_namespace)] * len(containing_tree.taxon_namespace)
    pvals = [1.0/len(containing_tree.taxon_namespace)] * len(containing_tree.taxon_namespace)
    sample_counts = multinomial_sample(n, pvals)
    leaf_nodes = containing_tree.leaf_nodes()
    assert len(sample_counts) == len(leaf_nodes)
    for sample_count, leaf_node in zip(sample_counts, leaf_nodes):
        leaf_node.num_individuals_sampled = sample_count
    filter_fn = lambda nd: nd.num_individuals_sampled > 0
    taxa_to_drop = set()
    dropped_nodes = containing_tree.filter_leaf_nodes(
            filter_fn,
            recursive=False)
    for nd in dropped_nodes:
        containing_tree.taxon_namespace.remove_taxon(nd.taxon)
    dropped_nodes = containing_tree.prune_leaves_without_taxa(recursive=True)
    return containing_tree

def generate_contained_trees(
        containing_tree,
        contained_taxon_namespace=None,
        population_size=1,
        total_number_of_individuals=200,
        num_gene_trees=5,
        rng=None):
    if contained_taxon_namespace is None:
        contained_taxon_namespace = dendropy.TaxonNamespace()
    contained_to_containing_map = {}
    assert len(containing_tree.taxon_namespace) > 0
    containing_tree = process_containing_tree_for_gene_samples(
            containing_tree=containing_tree,
            total_number_of_individuals=total_number_of_individuals,
            rng=rng)
    containing_tree_leaf_nodes = containing_tree.leaf_nodes()
    for sp_idx, sp_node in enumerate(containing_tree_leaf_nodes):
        sp_tax = sp_node.taxon
        for gidx in range(sp_node.num_individuals_sampled):
            glabel = "{sp}_{ind}^{sp}".format(sp=sp_tax.label, ind=gidx+1)
            # glabel = "{sp}^{sp}_{ind}".format(sp=sp_tax.label, ind=gidx+1)
            g = contained_taxon_namespace.require_taxon(label=glabel)
            g.population_label = sp_tax.label
            contained_to_containing_map[g] = sp_tax
    ct = reconcile.ContainingTree(
            containing_tree=containing_tree,
            contained_taxon_namespace=contained_taxon_namespace,
            contained_to_containing_taxon_map=contained_to_containing_map)
    gene_trees = dendropy.TreeList(taxon_namespace=contained_taxon_namespace)
    for gtidx in range(num_gene_trees):
        gt = ct.embed_contained_kingman(
                default_pop_size=population_size,
                rng=rng)
        gene_trees.append(gt)
    return containing_tree, gene_trees


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("source_trees",
            metavar="SOURCE_TREEFILE [SOURCE_TREEFILE [SOURCE_TREEFILE]]",
            nargs="+",
            help="Path to source of tree files. Specify '-' to read from standard input.")
    parser.add_argument("-f", "--input-format",
            default="nexus",
            dest="schema",
            help="Input trees format (default: $(default)s).")
    parser.add_argument("-z", "--random-seed",
            type=int,
            default=None,
            help="Seed for random number generator engine.")
    parser.add_argument("-t", "--title",
            default="bpprun",
            help="Run title (default: '%(default)s')")
    data_options = parser.add_argument_group("Data Options")
    data_options.add_argument("--population-size",
            type=int,
            default=10000,
            help="Population size (default: %(default)s).")
    data_options.add_argument("--total-number-of-individuals",
            type=int,
            default=200,
            help="Number of individuals sampled across all populations (default: %(default)s).")
    data_options.add_argument("--num-loci-per-individual",
            type=int,
            default=10,
            help="Number of loci sampled per individual (default: %(default)s).")
    data_options.add_argument("--num-characters-per-locus",
            type=int,
            default=1000,
            help="Number of characters sampled per locus (default: %(default)s).")
    data_options.add_argument("--mutation-rate-per-site",
            type=float,
            default=0.00001,
            help="Per-site mutation rate (default: %(default)s).")
    args = parser.parse_args()

    if args.random_seed is None:
        random_seed = random.randint(0, sys.maxint-1)
    else:
        random_seed = args.random_seed
    rng = random.Random(random_seed)
    s00.log("Random seed: {}".format(random_seed))

    sg = seqgen.SeqGen()
    sg.scale_branch_lengths = args.mutation_rate_per_site

    if "-" in args.source_trees:
        filepaths = sys.stdin.read().split("\n")
        args.source_trees.remove("-")
    else:
        filepaths = []

    manifest_entries = []
    filepaths.extend(args.source_trees)
    for idx, filepath in enumerate(filepaths):
        job_title = "{}_{:05d}".format(args.title, idx+1)
        manifest_entry = collections.OrderedDict()
        s00.log("{} of {}: {}: {}".format(idx+1, len(filepaths), job_title, filepath))
        try:
            source_tree = dendropy.Tree.get(
                    path=filepath,
                    schema=args.schema,
                    extract_comment_metadata=True,
                    preserve_underscores=True,
                    )
        except OSError, dendropy.DataError:
            s00.log("Skipping failed file: {}".format(filepath))
            continue

        manifest_entry["speciation_initiation_from_orthospecies_rate"] = try_to_coerce_to_float(source_tree.annotations["speciation_initiation_from_orthospecies_rate"].value)
        manifest_entry["speciation_initiation_from_incipient_species_rate"] = try_to_coerce_to_float(source_tree.annotations["speciation_initiation_from_incipient_species_rate"].value)
        manifest_entry["speciation_completion_rate"] = try_to_coerce_to_float(source_tree.annotations["speciation_completion_rate"].value)
        manifest_entry["orthospecies_extinction_rate"] = try_to_coerce_to_float(source_tree.annotations["orthospecies_extinction_rate"].value)
        manifest_entry["incipient_species_extinction_rate"] = try_to_coerce_to_float(source_tree.annotations["incipient_species_extinction_rate"].value)
        manifest_entry["max_time"] = try_to_coerce_to_float(source_tree.annotations["max_time"].value)
        manifest_entry["max_extant_orthospecies"] = try_to_coerce_to_float(source_tree.annotations["max_extant_orthospecies"].value)
        manifest_entry["num_extant_lineages"] = try_to_coerce_to_float(source_tree.annotations["num_extant_lineages"].value)
        manifest_entry["num_extant_orthospecies"] = try_to_coerce_to_float(source_tree.annotations["num_extant_orthospecies"].value)
        manifest_entry["source_tree_type"] = source_tree.annotations["tree_type"].value
        manifest_entry["population_size"] = args.population_size
        manifest_entry["total_number_of_individuals"] = args.total_number_of_individuals
        manifest_entry["num_loci_per_individual"] = args.num_loci_per_individual
        manifest_entry["mutation_rate_per_site"] = args.mutation_rate_per_site

        source_tree.calc_node_ages()
        original_containing_tree_num_species = len(source_tree.taxon_namespace)
        original_containing_tree_species_labels = " ".join(t.label for t in source_tree.taxon_namespace)

        containing_tree, gene_trees = generate_contained_trees(
                containing_tree=source_tree,
                total_number_of_individuals=args.total_number_of_individuals,
                num_gene_trees=args.num_loci_per_individual,
                population_size=args.population_size,
                rng=rng,
                )

        imap_filepath = "{}.input.imap.txt".format(job_title)
        f = open(imap_filepath, "w")
        for taxon in gene_trees.taxon_namespace:
            f.write("{}    {}\n".format(taxon.label.split("^")[1], taxon.population_label))
        f.write("\n//end of file")

        d0 = sg.generate(gene_trees)
        chars_filepath = "{}.input.chars.txt".format(job_title)
        f = open(chars_filepath, "w")
        for cm in d0.char_matrices:
            d0.write(file=f, schema="phylip")
            f.write("\n")

        out_filepath = "{}.results.out.txt".format(job_title)
        mcmc_filepath = "{}.results.mcmc.txt".format(job_title)
        final_containing_tree_num_species = len(containing_tree.taxon_namespace)

        # final_containing_tree_num_species = len(containing_tree.taxon_namespace)
        # final_containing_tree_species_labels = " ".join(t.label for t in containing_tree.taxon_namespace)
        # num_individuals_per_species = " ".join(str(args.num_individuals_per_population) for i in range(len(source_tree.taxon_namespace)))
        final_containing_tree_species_labels = []
        final_containing_tree_num_individuals_per_species = []
        for nd in containing_tree.leaf_node_iter():
            final_containing_tree_species_labels.append(nd.taxon.label)
            final_containing_tree_num_individuals_per_species.append(nd.num_individuals_sampled)
        final_containing_tree_num_species = len(final_containing_tree_species_labels)
        final_containing_tree_species_labels = " ".join(final_containing_tree_species_labels)
        final_containing_tree_num_individuals_per_species = " ".join([str(i) for i in final_containing_tree_num_individuals_per_species])

        theta_prior_mean = args.population_size * 4 * args.mutation_rate_per_site
        theta_prior_a = 2.0
        theta_prior_b = theta_prior_a/theta_prior_mean
        tau_prior_mean = containing_tree.seed_node.age
        tau_prior_a = 2.0
        tau_prior_b = tau_prior_a/tau_prior_mean

        manifest_entry["num_input_lineages"] = final_containing_tree_num_species
        manifest_entry["theta"] = theta_prior_mean
        manifest_entry["theta_prior_a"] = theta_prior_a
        manifest_entry["theta_prior_b"] = theta_prior_b
        manifest_entry["root_age"] = containing_tree.seed_node.age
        manifest_entry["tau_prior_a"] = tau_prior_a
        manifest_entry["tau_prior_b"] = tau_prior_b

        species_tree = containing_tree.as_string(
                schema="newick",
                suppress_leaf_taxon_labels=False,
                suppress_leaf_node_labels=True,
                suppress_internal_taxon_labels=True,
                suppress_internal_node_labels=True,
                suppress_rooting=True,
                suppress_edge_lengths=True,
                unquoted_underscores=True,
                preserve_spaces=True,
                store_tree_weights=False,
                suppress_annotations=True,
                suppress_item_comments=True,
                )
        bpp_config = BPP_TEMPLATE.format(
                chars_filepath=chars_filepath,
                imap_filepath=imap_filepath,
                out_filepath=out_filepath,
                mcmc_filepath=mcmc_filepath,
                num_species=final_containing_tree_num_species,
                species_labels=final_containing_tree_species_labels,
                num_individuals_per_species=final_containing_tree_num_individuals_per_species,
                species_tree=species_tree,
                theta_prior_a=theta_prior_a,
                theta_prior_b=theta_prior_b,
                tau_prior_a=tau_prior_a,
                tau_prior_b=tau_prior_b,
                num_loci=args.num_loci_per_individual,
                )
        bpp_ctl_filepath = "{}.input.bpp.ctl".format(job_title)
        f = open(bpp_ctl_filepath, "w")
        f.write(bpp_config)
        f.write("\n")

        jobf = open("{}.job.sge".format(job_title), "w")
        jobf.write("#! /bin/bash\n")
        jobf.write("#$ -cwd\n")
        jobf.write("#$ -V\n")
        jobf.write("#$ -S /bin/bash\n")
        jobf.write("#$ -l h_vmem=12G\n")
        jobf.write("#$ -l virtual_free=12G\n")
        jobf.write("bpp {}\n".format(bpp_ctl_filepath))

        manifest_entry["source_tree_path"] = filepath
        manifest_entry["results_filepath"] = out_filepath
        manifest_entry["mcmc_filepath"] = mcmc_filepath
        manifest_entries.append(manifest_entry)

    out = s00.open_output_file_for_csv_writer(
            filepath="{}_manifest.csv".format(args.title),
            append=False)
    with out:
        writer = csv.DictWriter(
                out,
                fieldnames=manifest_entries[0].keys(),
                restval="NA",
                delimiter=",",
                lineterminator=os.linesep,
                )
        writer.writeheader()
        writer.writerows(manifest_entries)

if __name__ == "__main__":
    main()


