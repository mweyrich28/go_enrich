package org.go;

import org.go.util.FileUtils;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.util.ArrayList;
import java.util.HashMap;

public class DAG
{

    private static final Logger logger = LoggerFactory.getLogger(DAG.class);
    private HashMap<String, GOEntry> nodeMap;
    private String type;
    private GOEntry root = null;


    public DAG(String type)
    {
        this.type = type;
        this.nodeMap = new HashMap<>();
    }

    public void insertNode(GOEntry go)
    {
        if (!nodeMap.containsKey(go.getId()))
        {
            nodeMap.put(go.getId(), go);
        }
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
