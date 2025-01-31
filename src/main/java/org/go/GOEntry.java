package org.go;

import java.util.*;

public class GOEntry
{

    private String id;
    private String annot;
    private boolean isTrue;
    private ArrayList<GOEntry> parents;
    private ArrayList<GOEntry> children;
    private HashSet<String> geneSymbols;
    private HashSet<String> reachableGoIDs;
    private HashMap<GOEntry, LinkerClass> highwayMap;
    private int depth;

    public GOEntry(String id, String annot)
    {
        this.id = id;
        this.annot = annot;
        this.geneSymbols = new HashSet<>();
        this.reachableGoIDs = new HashSet<>();
        this.reachableGoIDs.add(id);
        this.parents = new ArrayList<>();
        this.children = new ArrayList<>();
        this.highwayMap = new HashMap<>();
    }

    public void setAnnot(String annot)
    {
        if (this.annot == null)
        {
            this.annot = annot;
        }
    }

    public String getId()
    {
        return id;
    }

    public void addSymbol(String geneId)
    {
        this.geneSymbols.add(geneId);
    }

    public void addParent(GOEntry parent)
    {
        this.parents.add(parent);
    }

    public void addChild(GOEntry child)
    {
        this.children.add(child);
    }

    public ArrayList<GOEntry> getChildren()
    {
        return children;
    }

    @Override
    public String toString()
    {
        return this.depth + "\t" + this.id + "\t" + this.annot + "\t" + children.size();
    }

    public void setTrue()
    {
        isTrue = true;
    }

    public void inheritGeneSymbols()
    {
        for (GOEntry child : children)
        {
            child.inheritGeneSymbols();
            this.getGeneSymbols().addAll(child.geneSymbols);
        }
    }

    public void inheritGoIds()
    {
        for (GOEntry child : children)
        {
            child.inheritGoIds();
            this.reachableGoIDs.addAll(child.reachableGoIDs);
        }
    }

    public HashSet<String> getGeneSymbols()
    {
        return geneSymbols;
    }

    public void setDepth(int depth)
    {
        this.depth = depth;
    }

    public void passDepth()
    {
        int newDepth = this.depth;
        newDepth++;
        for (GOEntry child : children)
        {
            child.setDepth(newDepth);
            child.passDepth();
        }
    }

    public String getAnnot()
    {
        return annot;
    }

    public boolean isTrue()
    {
        return isTrue;
    }

    public String getShortestPathToTrue(HashSet<String> groundTruth, DAG dag)
    {
        if (groundTruth == null || groundTruth.isEmpty() || groundTruth.contains(this.id))
        {
            return "";
        }
        GOEntry target = null;
        int c = Integer.MAX_VALUE;
        for (String trueGOID : groundTruth)
        {
            GOEntry trueGO = dag.getNodeMap().get(trueGOID);
            int length = this.highwayMap.get(trueGO).getK();
            if (length < c)
            {
                c = length;
                target = trueGO;
            }
        }


        GOEntry curr = this;
        StringBuilder sb = new StringBuilder();
        boolean seenLCA = false;
        while (curr.highwayMap.containsKey(target))
        {
            if (curr.equals(this))
            {
                sb.append(curr.annot);
            }
            // LCA
            else if (curr.reachableGoIDs.contains(target.id) && !seenLCA)
            {
                seenLCA = true;
                sb.append("|").append(curr.annot).append(" * ");
                if (curr.id.equals(target.id))
                {
                    break;
                }
            }
            else
            {
                sb.append("|").append(curr.annot);
            }

            curr = curr.highwayMap.get(target).getGo();
            if (curr.id.equals(target.id))
            {
                if (seenLCA)
                {
                    sb.append("|").append(curr.annot);
                }
                else
                {
                    sb.append("|").append(curr.annot).append(" * ");
                }
                break;
            }

            // extra check
            if (groundTruth.contains(curr.id))
            {
                sb.append("|").append(curr.annot).append(" * ");
                break;
            }
        }

        return sb.toString();
    }

    public void signalShortestPathUp(int k, GOEntry trueGo)
    {
        k++;
        for (GOEntry parent : parents)
        {
            if (!parent.highwayMap.containsKey(trueGo))
            {
                LinkerClass link = new LinkerClass(this, k);
                parent.highwayMap.put(trueGo, link);
            }
            else
            {
                if (parent.highwayMap.get(trueGo).getK() > k)
                {
                    parent.highwayMap.get(trueGo).setGo(this);
                    parent.highwayMap.get(trueGo).setK(k);
                }
            }
            parent.signalShortestPathUp(k, trueGo);
        }
    }

    public void signalShortestPathDown(int k, GOEntry trueGo)
    {
        k++;
        for (GOEntry child : this.children)
        {
            if (!child.highwayMap.containsKey(trueGo))
            {
                LinkerClass link = new LinkerClass(this, k);
                child.highwayMap.put(trueGo, link);
            }
            else
            {
                if (child.highwayMap.get(trueGo).getK() > k)
                {
                    child.highwayMap.get(trueGo).setGo(this);
                    child.highwayMap.get(trueGo).setK(k);
                }
            }
            child.signalShortestPathDown(k, trueGo);
        }
    }
    public void propagateShortestPaths(GOEntry trueGo)
    {
        LinkerClass currentPath = highwayMap.get(trueGo);
        if (currentPath != null)
        {
            GOEntry nextNode = currentPath.getGo();
            if (nextNode.highwayMap.containsKey(trueGo))
            {
                int pathThroughLinker = nextNode.highwayMap.get(trueGo).getK() + 1;
                if (pathThroughLinker < currentPath.getK())
                {
                    currentPath.setK(pathThroughLinker);
                }
            }
        }

        int currentK = currentPath != null ? currentPath.getK() : Integer.MAX_VALUE;
        for (GOEntry child : children)
        {
            LinkerClass childPath = child.highwayMap.get(trueGo);
            int childK = childPath != null ? childPath.getK() : Integer.MAX_VALUE;

            if (currentK != Integer.MAX_VALUE && currentK + 1 < childK)
            {
                if (childPath == null)
                {
                    childPath = new LinkerClass(this, currentK + 1);
                    child.highwayMap.put(trueGo, childPath);
                }
                else
                {
                    childPath.setGo(this);
                    childPath.setK(currentK + 1);
                }
            }
            child.propagateShortestPaths(trueGo);
        }
    }

    public HashSet<String> getReachableGoIDs()
    {
        return reachableGoIDs;
    }

    public ArrayList<GOEntry> getParents()
    {
        return parents;
    }

}
