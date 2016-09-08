#! /usr/bin/env python

import os
import sys
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, SCRIPT_DIR)
import s00

import random
import argparse
from dendropy.model import protractedspeciation
from dendropy.model import reconcile
import dendropy

def describe_tree(out, tree, title):
    if title:
        out.write("{}:\n".format(title))
    out.write("   -  ")
    tree.write(file=out, schema="newick", suppress_annotations=False)
    out.write("   -  Number of tips: {}\n".format(len(tree.leaf_nodes())))
    out.write("   -        Root age: {}\n".format(tree.seed_node.age))

def set_lineage_tree_taxa_from_labels(tree):
    seen = set()
    for idx, nd in enumerate(tree):
        # if nd.annotations["is_full_speciation_event"]:
        #     prefix = "S"
        # else:
        #     prefix = "I"
        # prefix = "Q"
        # label = "{}{:02d}".format(prefix, idx+1)
        if nd.is_leaf():
            label = nd.label
            assert label not in seen
            seen.add(label)
            nd.taxon = tree.taxon_namespace.require_taxon(label=label)
        else:
            nd.label = None
        # else:
        #     nd.label = label

def set_orthospecies_tree_taxa_from_labels(tree):
    seen = set()
    for idx, nd in enumerate(tree):
        # prefix = "S"
        # label = "{}{:02d}".format(prefix, idx+1)
        # if nd.is_leaf():
        #     if not nd.included_lineage_tree_leaf_nodes:
        #         print("Not found: {}".format(id(nd)))
        #     lineage_label_set = [i.label for i in nd.included_lineage_tree_leaf_nodes]
        #     label = "_".join(lineage_label_set)
        #     assert label
        # else:
        #     label = "{}{:02d}".format(prefix, idx+1)
        # if label in seen:
        #     print("{}: {}".format(label, seen))
        if nd.is_leaf():
            label = nd.label
            assert label not in seen, label
            seen.add(label)
            nd.taxon = tree.taxon_namespace.require_taxon(label=label)
        else:
            nd.label = None
        # nd.annotations["included_lineages"] = [i.label for i in nd.included_lineage_tree_leaf_nodes]
        # nd.annotations["num_included_lineages"] = len(nd.included_lineage_tree_leaf_nodes)

def main():
    parser = argparse.ArgumentParser()
    parameter_options = parser.add_argument_group("Model Parameters")
    parameter_options.add_argument("--b1", "--speciation_initiation_from_orthospecies_rate",
            type=float,
            dest="speciation_initiation_from_orthospecies_rate",
            default=0.5,
            help="Rate at which orthospecies give rise to new incipient species [default: %(default)s].")
    parameter_options.add_argument("--b2", "--speciation_initiation_from_incipient_species_rate",
            type=float,
            dest="speciation_initiation_from_incipient_species_rate",
            default=0.5,
            help="Rate at which incipient species give rise to new incipient species [default: %(default)s].")
    parameter_options.add_argument("--c1", "--speciation-completion-rate",
            type=float,
            dest="speciation_completion_rate",
            default=0.1,
            help="Rate at which incipient species become orthospecies [default: %(default)s].")
    parameter_options.add_argument("--e1", "--orthospecies-extinction-rate",
            type=float,
            dest="orthospecies_extinction_rate",
            default=0.1,
            help="Rate at which orthospecies go extinct [default: %(default)s].")
    parameter_options.add_argument("--e2", "--incipient-species-extinction-rate",
            type=float,
            dest="incipient_species_extinction_rate",
            default=0.1,
            help="Rate at which incipient species go extinct [default: %(default)s].")
    termination_options = parser.add_argument_group("Simulation Termination Conditions")
    termination_options.add_argument("--max-time",
            type=float,
            default=None,
            help="Maximum length of time to to run (default: %(default)s).")
    termination_options.add_argument("--max-extant-orthospecies",
            type=int,
            default=None,
            help="Maximum number of orthospecies to generate (default: %(default)s).")
    termination_options.add_argument("--max-extant-lineages",
            type=int,
            default=None,
            help="Maximum number of lineages to generate (default: %(default)s).")
    # data_options = parser.add_argument_group("Data Options")
    # data_options.add_argument("--population-size",
    #         type=int,
    #         default=10000,
    #         help="Population size (default: %(default)s).")
    # data_options.add_argument("--num-individuals-per-population",
    #         type=int,
    #         default=4,
    #         help="Number of individuals sampled per incipient species lineage (default: %(default)s).")
    # data_options.add_argument("--num-loci-per-individual",
    #         type=int,
    #         default=10,
    #         help="Number of loci sampled per individual (default: %(default)s).")
    # data_options.add_argument("--num-characters-per-locus",
    #         type=int,
    #         default=10,
    #         help="Number of characters sampled per locus (default: %(default)s).")
    # data_options.add_argument("--mutation-rate-per-site",
    #         type=float,
    #         default=0.00001,
    #         help="Per-site mutation rate (default: %(default)s).")
    run_options = parser.add_argument_group("Run Options")
    run_options.add_argument("-t", "--title",
            default="psmrun",
            help="Run title (default: '%(default)s')")
    run_options.add_argument("-n", "--nreps",
            type=int,
            default=10,
            help="Number of replicates (default: %(default)s).")
    run_options.add_argument("-z", "--random-seed",
            type=int,
            default=None,
            help="Seed for random number generator engine.")
    args = parser.parse_args()

    if not args.max_time and not args.max_extant_orthospecies and not args.max_extant_lineages:
        sys.exit("Need to specify termination condition, at least one of: '--max-time', '--max-extant-orthospecies', '--max-extant-lineages'")
    if args.random_seed is None:
        random_seed = random.randint(0, sys.maxint-1)
    else:
        random_seed = args.random_seed

    s00.log("Random seed: {}".format(random_seed))
    rng = random.Random(random_seed)
    psm = protractedspeciation.ProtractedSpeciationProcess(
            speciation_initiation_from_orthospecies_rate=args.speciation_initiation_from_orthospecies_rate,
            speciation_initiation_from_incipient_species_rate=args.speciation_initiation_from_incipient_species_rate,
            speciation_completion_rate=args.speciation_completion_rate,
            orthospecies_extinction_rate=args.orthospecies_extinction_rate,
            incipient_species_extinction_rate=args.incipient_species_extinction_rate,
            rng=rng,)
    # sg = seqgen.SeqGen()
    # sg.scale_branch_lengths = args.mutation_rate_per_site
    for rep in range(args.nreps):
        job_title = "{}_{:05d}".format(args.title, rep+1)
        s00.log("Replicate {} of {}: {}".format(rep+1, args.nreps, job_title))
        while True:
            lineage_tree, orthospecies_tree = psm.generate_sample(
                    max_time=args.max_time,
                    max_extant_orthospecies=args.max_extant_orthospecies,
                    max_extant_lineages=args.max_extant_lineages,
                    is_initial_lineage_orthospecies=True,
                    # is_correlate_lineage_and_species_trees=True,
                    )
            num_extant_lineages = len([nd for nd in lineage_tree.leaf_node_iter()])
            num_extant_orthospecies = len([nd for nd in orthospecies_tree.leaf_node_iter()])
            if num_extant_lineages < 3:
                s00.log("Too few lineages ({}): Repeating replicate {} of {}: {}".format(num_extant_lineages, rep+1, args.nreps, job_title))
            # elif num_extant_orthospecies < 2:
            #     s00.log("Too few orthospecies ({}): Repeating replicate {} of {}: {}".format(num_extant_orthospecies, rep+1, args.nreps, job_title))
            else:
                break
        set_lineage_tree_taxa_from_labels(lineage_tree)
        set_orthospecies_tree_taxa_from_labels(orthospecies_tree)
        s00.log("Lineage tree size: {}, Orthospecies tree size: {}".format(num_extant_lineages, num_extant_orthospecies))
        for tree in (lineage_tree, orthospecies_tree):
            tree.annotations["speciation_initiation_from_orthospecies_rate"] = args.speciation_initiation_from_orthospecies_rate
            tree.annotations["speciation_initiation_from_incipient_species_rate"] = args.speciation_initiation_from_incipient_species_rate
            tree.annotations["speciation_completion_rate"] = args.speciation_completion_rate
            tree.annotations["orthospecies_extinction_rate"] = args.orthospecies_extinction_rate
            tree.annotations["incipient_species_extinction_rate"] = args.incipient_species_extinction_rate
            tree.annotations["max_simulation_time"] = args.max_time
            tree.annotations["max_extant_orthospecies"] = args.max_extant_orthospecies
            tree.annotations["num_extant_lineages"] = num_extant_lineages
            tree.annotations["num_extant_orthospecies"] = num_extant_orthospecies
        lineage_tree.annotations["tree_type"] = "lineage"
        orthospecies_tree.annotations["tree_type"] = "orthospecies"
        lineage_tree.write(
                path="{}.lineages.tre".format(job_title),
                schema="nexus",
                suppress_annotations=False)
        orthospecies_tree.write(
                path="{}.species.tre".format(job_title),
                schema="nexus",
                suppress_annotations=False)

def __STASH__():
    """
        # lineage_tree.calc_node_ages()
        # orthospecies_tree.calc_node_ages()
        # _log("    Incipient species tree: {} tips, root age = {} ({} mutation units)".format(len(lineage_tree.leaf_nodes()), lineage_tree.seed_node.age, lineage_tree.seed_node.age * args.mutation_rate_per_site,))
        # _log("    Orthospecies tree:      {} tips, root age = {} ({} mutation units)".format(len(orthospecies_tree.leaf_nodes()), orthospecies_tree.seed_node.age))
        set_lineage_tree_taxa_from_labels(lineage_tree)
        set_orthospecies_tree_taxa_from_labels(orthospecies_tree)
        lineage_tree.write(path="x1.nexus", schema="nexus")
        orthospecies_tree.write(path="x2.nexus", schema="nexus")

        logf = open("{}.setup.log".format(job_title), "w")
        logf.write("-  Replicate {} of {} generated by command: {}\n".format(rep+1, args.nreps, " ".join(sys.argv)))
        logf.write("\n")
        logf.write("-  Random seed used: {}\n".format(random_seed))
        logf.write("\n")
        describe_tree(logf, lineage_tree, "-  Lineage Tree Profile")
        logf.write("\n")
        describe_tree(logf, orthospecies_tree, "-  Orthospecies Tree Profile")
        logf.write("\n")
        logf.write("-  Protracted Speciation Model Parameters\n")
        logf.write("   -       Speciation initiation from orthospecies rate: {}\n".format(args.speciation_initiation_from_orthospecies_rate))
        logf.write("   -  Speciation initiation from incipient species rate: {}\n".format(args.speciation_initiation_from_incipient_species_rate))
        logf.write("   -                         Speciation completion rate: {}\n".format(args.speciation_completion_rate))
        logf.write("   -                       Orthospecies extinction rate: {}\n".format(args.orthospecies_extinction_rate))
        logf.write("   -                  Incipient species extinction rate: {}\n".format(args.incipient_species_extinction_rate))
        logf.write("   -               Termination: Maximum simulation time: {}\n".format(args.max_time))
        logf.write("   -                  Termination: Maximum orthospecies: {}\n".format(args.max_extant_orthospecies))
        logf.write("   -                      Termination: Maximum lineages: {}\n".format(args.max_extant_lineages))
        logf.write("\n")
        logf.write("-  Data Generation Parameters\n")
        logf.write("   -                                    Population size: {}\n".format(args.population_size))
        logf.write("   -                 Individuals per species/population: {}\n".format(args.num_individuals_per_population))
        logf.write("   -                      Number of loci per individual: {}\n".format(args.num_loci_per_individual))
        logf.write("   -                             Per-site mutation rate: {}\n".format(args.mutation_rate_per_site))
        logf.write("\n")
        # logf.write("-  Analysis Setup\n")
        # logf.write("   -                          Theta prior [a, b (mean)]: {}, {} ({})\n".format(theta_prior_a, theta_prior_b, theta_prior_mean))
        # logf.write("   -                            Tau prior [a, b (mean)]: {}, {} ({})".format(tau_prior_a, tau_prior_b, tau_prior_mean))
        # logf.write("\n")


        gene_trees = generate_contained_trees(
                containing_tree=lineage_tree,
                num_individuals_per_population=args.num_individuals_per_population,
                num_gene_trees=args.num_loci_per_individual,
                population_size=args.population_size,
                rng=rng,
                # taxon_namespace=gene_trees_taxon_namespace,
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
        num_species = len(lineage_tree.taxon_namespace)
        species_labels = " ".join(t.label for t in lineage_tree.taxon_namespace)
        num_individuals_per_species = " ".join(str(args.num_individuals_per_population) for i in range(len(lineage_tree.taxon_namespace)))
        theta_prior_mean = args.population_size * 4 * args.mutation_rate_per_site
        theta_prior_a = 2.0
        theta_prior_b = theta_prior_a/theta_prior_mean
        tau_prior_mean = lineage_tree.seed_node.age
        tau_prior_a = 2.0
        tau_prior_b = tau_prior_a/tau_prior_mean
        species_tree = lineage_tree.as_string(
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
                num_species=num_species,
                species_labels=species_labels,
                num_individuals_per_species=num_individuals_per_species,
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

        lineage_tree.write(path="{}.setup.lineage-tree.nexus".format(job_title), schema="nexus")
        orthospecies_tree.write(path="{}.setup.confirmed-species-tree.nexus".format(job_title), schema="nexus")
        """


if __name__ == "__main__":
    main()


