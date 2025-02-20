package org.go;

import net.sourceforge.argparse4j.ArgumentParsers;
import net.sourceforge.argparse4j.inf.ArgumentParser;
import net.sourceforge.argparse4j.inf.ArgumentParserException;
import net.sourceforge.argparse4j.inf.Namespace;

import java.io.IOException;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.logging.Logger;

import org.go.util.FileUtils;
import org.slf4j.LoggerFactory;

public class Main
{

    public static HashSet<String> VALIDTYPES = new HashSet<>();
    public static HashSet<String> VALIDROOTS = new HashSet<>();

    private static final org.slf4j.Logger logger = LoggerFactory.getLogger(Main.class);

    public static void main(String[] args) throws IOException
    {
        ArgumentParser parser = ArgumentParsers.newFor("GoEnrich").build().defaultHelp(true).description("");
        try
        {
            double start = System.currentTimeMillis();
            logger.info("Starting with GO Enrichment Assignment 5");
            parser.addArgument("-obo").required(true).help("Path to obo file.");
            parser.addArgument("-root").required(true).help("Type of Ontology (molecular_function|biological_process|cellular_component).");
            parser.addArgument("-mapping").required(true).help("Path to mapping file (SYMBOL -> GO_ID).");
            parser.addArgument("-mappingtype").required(true).help("Format of the mapping-File (go|ensembl).");
            parser.addArgument("-overlapout").help("Information about DAG entries with shared mapped genes is written into this file.");
            parser.addArgument("-enrich").required(true).help("Path to enrichment analysis file.");
            parser.addArgument("-o").required(true).help("Path to output file.");
            parser.addArgument("-minsize").required(true).help("Min amount of genes per GO entry.");
            parser.addArgument("-maxsize").required(true).help("Max amount of genes per GO entry.");

            Namespace ns = parser.parseArgs(args);

            String obo = ns.getString("obo");
            String root = ns.getString("root");
            String mapping = ns.getString("mapping");
            String mappingType = ns.getString("mappingtype");
            String overlapOut = ns.getString("overlapout");
            String enrich = ns.getString("enrich");
            String out = ns.getString("o");
            int min = Integer.parseInt(ns.getString("minsize"));
            int max = Integer.parseInt(ns.getString("maxsize"));

            // validate params
            validateParams(root, mappingType);


            logger.info("Parsing input files...");
            double s2 = System.currentTimeMillis();
            // init DAG
            DAG dag = FileUtils.parseOBO(obo, root);

            // read mapping
            if (mappingType.equals("go"))
            {
                FileUtils.parseGAF(mapping, dag);
            }
            else
            {
                FileUtils.parseENSEMBL(mapping, dag);
            }
            logger.info(String.format("Total time needed for parsing files: %s seconds", (System.currentTimeMillis() - s2) / 1000.0));

            // read enriched genes
            HashMap<String, Gene> enrichedGeneMap = FileUtils.parseEnrichedGenes(enrich, dag);

            // fully init DAG
            logger.info("Starting to optimize DAG...");
            s2 = System.currentTimeMillis();
            dag.propagateGenes();
            dag.propagateIds();
            dag.optimizePathsToTrue();
            logger.info(String.format("Total time needed for optimizing DAG: %s seconds", (System.currentTimeMillis() - s2) / 1000.0));

            // init main analysis class
            EnrichmentAnalysis enrichmentAnalysis = new EnrichmentAnalysis(dag, enrichedGeneMap, min, max, out);
            //            enrichmentAnalysis.analyze();
            enrichmentAnalysis.analyzeParallelized();


            // extra analysis for report
            //            dag.initLeafNodes();
            //            System.out.println("Num genes " + dag.getRoot().getGeneSymbols().size());
            //            System.out.println("Num gene sets " + dag.getNodeMap().size());
            //            System.out.println("Num leafs " + dag.computeNumLeafs());
            //            dag.computePathLengths();
            //            dag.writeGeneSetSizes(mappingType + "_goSizes.tsv");
            //            dag.writeDiffSize(mappingType + "_goDiff.tsv");

            if (overlapOut != null)
            {
                //                 enrichmentAnalysis.goFeatures(overlapOut);
                enrichmentAnalysis.goFeaturesParallelized(overlapOut);
            }
            logger.info(String.format("Total time needed for Assignemnt 5: %s seconds", (System.currentTimeMillis() - start) / 1000.0));
        }
        catch (ArgumentParserException e)
        {
            parser.printHelp();
        }
    }

    public static void validateParams(String root, String type)
    {
        VALIDROOTS.addAll(List.of(new String[]{"biological_process", "molecular_function", "cellular_component"}));
        VALIDTYPES.addAll(List.of(new String[]{"go", "ensembl"}));
        if (!VALIDTYPES.contains(type))
        {
            throw new RuntimeException("INVALID MAPPING TYPE! PROVIDE -mappingtype [go|ensembl]");
        }
        if (!VALIDROOTS.contains(root))
        {
            throw new RuntimeException("INVALID GO ONTOLOGY! PROVIDE -root [biological_process|molecular_function|cellular_component]");
        }
    }
}