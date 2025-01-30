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

    public String getShortestPathToTrueNewNew(HashSet<String> groundTruth, DAG dag)
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
                if (curr.id.equals(target.id)) {
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
                sb.append("|").append(curr.annot);
                break;
            }
        }
//        System.out.println(sb);

        //        GOEntry curr = LCA;
        //        while (visited.containsKey(curr))
        //        {
        //            if (curr.id.equals(LCA.id))
        //            {
        //                sb.append(curr.annot).append(" * ");
        //            }
        //            else
        //            {
        //                sb.insert(0, curr.annot + "|");
        //            }
        //            GOEntry p = visited.get(curr);
        //            curr = p;
        //        }
        //
        //        curr = LCA;
        //        while (curr.highwayMap.containsKey(dag.getNodeMap().get(closestTrueGO)))
        //        {
        //            curr = curr.highwayMap.get(dag.getNodeMap().get(closestTrueGO)).getGo();
        //            sb.append("|").append(curr.annot);
        //        }


        return sb.toString();
    }

    public String getShortestPathToTrueNew(HashSet<String> groundTruth, DAG dag)
    {
        if (groundTruth == null || groundTruth.isEmpty() || groundTruth.contains(this.id))
        {
            return "";
        }

        Queue<GOEntry> queue = new LinkedList<>();
        HashMap<GOEntry, GOEntry> visited = new HashMap<>();
        queue.add(this);
        visited.put(this, null);
        GOEntry LCA = null;
        ArrayList<String> potentialTruths = new ArrayList<>();
        String closestTrueGO = null;

        while (!queue.isEmpty() && LCA == null)
        {
            GOEntry current = queue.poll();
            for (GOEntry parent : current.parents)
            {
                if (!visited.containsKey(parent))
                {
                    visited.put(parent, current);
                    if (groundTruth.contains(parent.id))
                    {
                        LCA = parent;
                        break;
                    }

                    // check if ground trough can be reached
                    HashSet<String> copyTruth = new HashSet<>(groundTruth);
                    copyTruth.retainAll(parent.reachableGoIDs);

                    if (!copyTruth.isEmpty())
                    {
                        potentialTruths.addAll(copyTruth);
                        int pathLength = Integer.MAX_VALUE;
                        for (String potentialTruth : potentialTruths)
                        {
                            if (parent.highwayMap.get(dag.getNodeMap().get(potentialTruth)).getK() < pathLength)
                            {
                                LCA = parent;
                                pathLength = parent.highwayMap.get(dag.getNodeMap().get(potentialTruth)).getK();
                                closestTrueGO = potentialTruth;
                            }
                        }
                        break;
                    }
                    queue.add(parent);
                }
            }
        }

        if (LCA == null)
        {
            return "";
        }

        GOEntry curr = LCA;
        StringBuilder sb = new StringBuilder();
        while (visited.containsKey(curr))
        {
            if (curr.id.equals(LCA.id))
            {
                sb.append(curr.annot).append(" * ");
            }
            else
            {
                sb.insert(0, curr.annot + "|");
            }
            GOEntry p = visited.get(curr);
            curr = p;
        }

        curr = LCA;
        while (curr.highwayMap.containsKey(dag.getNodeMap().get(closestTrueGO)))
        {
            curr = curr.highwayMap.get(dag.getNodeMap().get(closestTrueGO)).getGo();
            sb.append("|").append(curr.annot);
        }


        return sb.toString();
    }

    public String getShortestPathToTrue(HashSet<String> groundTruth, DAG dag)
    {
        if (groundTruth == null || groundTruth.isEmpty() || groundTruth.contains(this.id))
        {
            return "";
        }
        //        if (this.id.equals("GO:1900371"))
        //        {
        //            System.out.println();
        //        }

        Queue<GOEntry> queue = new LinkedList<>();
        HashMap<GOEntry, GOEntry> visited = new HashMap<>();
        queue.add(this);
        visited.put(this, null);
        GOEntry LCA = null;
        ArrayList<String> potentialTruths = new ArrayList<>();

        while (!queue.isEmpty() && LCA == null)
        {
            GOEntry current = queue.poll();
            for (GOEntry parent : current.parents)
            {
                if (!visited.containsKey(parent))
                {
                    visited.put(parent, current);
                    if (groundTruth.contains(parent.id))
                    {
                        LCA = parent;
                        break;
                    }

                    // check if ground trough can be reached
                    HashSet<String> copyTruth = new HashSet<>(groundTruth);
                    copyTruth.retainAll(parent.reachableGoIDs);

                    if (!copyTruth.isEmpty())
                    {
                        potentialTruths.addAll(copyTruth);
                        LCA = parent;
                        break;
                    }
                    queue.add(parent);
                }
            }
        }

        if (LCA == null)
        {
            return "";
        }

        GOEntry curr = LCA;
        StringBuilder sb = new StringBuilder();
        while (visited.containsKey(curr))
        {
            if (curr.id.equals(LCA.id))
            {
                sb.append(curr.annot).append(" * ");
            }
            else
            {
                sb.insert(0, curr.annot + "|");
            }
            GOEntry p = visited.get(curr);
            curr = p;
        }

        String shortestPath = null;
        int shortestLength = Integer.MAX_VALUE;
        for (String potTruthId : potentialTruths)
        {
            GOEntry potGO = dag.getNodeMap().get(potTruthId);
            String pathCandidate = LCA.findPathToTruth(potGO);

            if (shortestPath == null)
            {
                shortestPath = pathCandidate;
                shortestLength = calcPathLength(pathCandidate);
                continue;
            }
            int pathLength = calcPathLength(pathCandidate);
            if (pathLength <= shortestLength)
            {
                shortestPath = pathCandidate;
                shortestLength = pathLength;
            }
        }
        if (shortestPath != null)
        {
            sb.append(shortestPath);
        }
        //        System.out.println(this.id);
        //        System.out.println(sb);


        return sb.toString();
    }

    public int calcPathLength(String path)
    {
        int c = 0;
        for (int i = 0; i < path.length(); i++)
        {
            if (path.charAt(i) == '|')
            {
                c++;
            }
        }
        return c;
    }


    public String findPathToTruth(GOEntry truth)
    {

        Queue<GOEntry> queue = new LinkedList<>();
        HashMap<GOEntry, GOEntry> visited = new HashMap<>();
        queue.add(this);
        visited.put(this, null);

        while (!queue.isEmpty())
        {
            GOEntry current = queue.poll();
            for (GOEntry child : current.children)
            {
                if (!visited.containsKey(child))
                {
                    visited.put(child, current);
                    if (truth.id.equals(child.id))
                    {

                        StringBuilder sb = new StringBuilder();
                        GOEntry curr = truth;
                        while (visited.containsKey(curr))
                        {
                            if (!curr.id.equals(this.id))
                            {
                                sb.append("|" + curr.annot);
                            }
                            GOEntry p = visited.get(curr);
                            curr = p;
                        }
                        return sb.toString();

                    }
                    queue.add(child);
                }
            }
        }
        return "";
    }

    public void signalShortestPath(int k, GOEntry trueGo)
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
            parent.signalShortestPath(k, trueGo);
        }
    }

    public HashMap<GOEntry, LinkerClass> getHighwayMap()
    {
        return highwayMap;
    }

    public void propagateShortestPaths(GOEntry trueGo)
    {
        // First, check if we have a path to this trueGo
        LinkerClass currentPath = highwayMap.get(trueGo);
        int currentK = currentPath != null ? currentPath.getK() : Integer.MAX_VALUE;

        // Check all children
        for (GOEntry child : children)
        {  // Assuming you have a children collection
            LinkerClass childPath = child.highwayMap.get(trueGo);
            int childK = childPath != null ? childPath.getK() : Integer.MAX_VALUE;

            // If we have a better path than the child's current path
            if (currentK != Integer.MAX_VALUE && currentK + 1 < childK)
            {
                // Either create new or update existing path
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

            // Always propagate to children to ensure complete coverage
            child.propagateShortestPaths(trueGo);
        }
    }
}
