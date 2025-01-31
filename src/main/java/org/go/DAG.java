package org.go;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;

public class DAG
{

    private static final Logger logger = LoggerFactory.getLogger(DAG.class);
    private HashMap<String, GOEntry> nodeMap;
    private String type;
    private GOEntry root = null;
    private HashSet<GOEntry> trueGoEntries;
    private HashSet<String> trueGoIds;


    public DAG(String type)
    {
        this.type = type;
        this.nodeMap = new HashMap<>();
        this.trueGoEntries = new HashSet<>();
        this.trueGoIds = new HashSet<>();
    }

    public void insertNode(GOEntry go)
    {
        if (!nodeMap.containsKey(go.getId()))
        {
            nodeMap.put(go.getId(), go);
        }
    }

    public void addTrueGoEntry(GOEntry go)
    {
        this.trueGoEntries.add(go);
    }

    public void addTrueGoId(String id)
    {
        this.trueGoIds.add(id);
    }

    public HashSet<String> getTrueGoIds()
    {
        return trueGoIds;
    }

    public HashSet<GOEntry> getTrueGoEntries()
    {
        return trueGoEntries;
    }

    @Override
    public String toString()
    {
        StringBuilder sb = new StringBuilder();
        for (String goId : nodeMap.keySet())
        {
            sb.append(goId);
            if (nodeMap.get(goId).getChildren().size() == 0)
            {
                sb.append("\t STOP");
            }
            else
            {
                sb.append("\tCHILDREN ->");
                sb.append("\t");
                ArrayList<String> ids = new ArrayList<>();
                for (GOEntry go : nodeMap.get(goId).getChildren())
                {
                    ids.add(go.getId());
                }
                sb.append(String.join(" | ", ids));
            }
            sb.append("\n\n");
        }

        return sb.toString();
    }

    public void setRoot(GOEntry root)
    {
        if (this.root == null)
        {
            this.root = root;
        }
    }


    public HashMap<String, GOEntry> getNodeMap()
    {
        return nodeMap;
    }

    public void propagateGenes()
    {
        logger.info("Propagating Symbols through DAG...");
        double start = System.currentTimeMillis();
        root.inheritGeneSymbols();
        logger.info(String.format("Time needed for propagating Symbols: %s seconds", (System.currentTimeMillis() - start) / 1000.0));
    }

    public void propagateIds()
    {
        logger.info("Propagating GO Ids through DAG...");
        double start = System.currentTimeMillis();
        root.inheritGoIds();
        logger.info(String.format("Time needed for propagating GO Ids: %s seconds", (System.currentTimeMillis() - start) / 1000.0));
    }

    public void optimizePathsToTrue()
    {
        logger.info("Optimizing paths to true GOs in DAG...");
        double start = System.currentTimeMillis();
        for(GOEntry truth: this.getTrueGoEntries())
        {
            truth.signalShortestPathUp(0, truth);
        }
        for(GOEntry truth: this.getTrueGoEntries())
        {
            root.propagateShortestPaths(truth);
        }
        for(GOEntry truth: this.getTrueGoEntries())
        {
            truth.signalShortestPathDown(0, truth);
        }

//
//        for(GOEntry truth: this.getTrueGoEntries())
//        {
//            root.checkDirectPaths(truth);
//        }
//
//        for(GOEntry truth: this.getTrueGoEntries())
//        {
//            root.updatePathDistances(truth);
//        }
        logger.info(String.format("Time needed for optimizing paths: %s seconds", (System.currentTimeMillis() - start) / 1000.0));
    }


    public void calculateDepth()
    {
        logger.info("Calculating depth...");
        double start = System.currentTimeMillis();
        this.root.setDepth(0);
        this.root.passDepth();
        logger.info(String.format("Time needed for calculating depth: %s seconds", (System.currentTimeMillis() - start) / 1000.0));
    }

    public GOEntry getRoot()
    {
        return root;
    }
}
