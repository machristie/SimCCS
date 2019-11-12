/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package solver;

import dataStore.DataInOut;
import dataStore.DataStorer;
import dataStore.Edge;
import dataStore.LinearComponent;
import dataStore.Sink;
import dataStore.Solution;
import dataStore.Source;
import dataStore.UnidirEdge;
import ilog.concert.*;
import ilog.cplex.*;
import java.io.File;
import java.text.DateFormat;
import java.text.SimpleDateFormat;
import java.util.Date;
import java.util.HashMap;
import java.util.HashSet;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 *
 * @author yaw
 */
public class FlowHeuristic {

    public static void run(DataStorer data, double crf, double numYears, double captureTarget, String basePath, String dataset, String scenario) {
        // create directory
        HashMap<Edge, Double> edgeHostingAmounts = new HashMap<>();
        DateFormat dateFormat = new SimpleDateFormat("ddMMyyy-HHmmssss");
        Date date = new Date();
        String run = "flowCap" + dateFormat.format(date);
        File solutionDirectory = new File(basePath + "/" + dataset + "/Scenarios/" + scenario + "/Results/" + run + "/");
        solutionDirectory.mkdir();

        int numIterations = 15;
        for (int i = 0; i < numIterations; i++) {
            //if (i % 10 == 0 || i == numIterations - 1) {
                System.out.print("Iteration " + i + ": ");
                edgeHostingAmounts = makeSolveLP(data, crf, numYears, captureTarget, edgeHostingAmounts, solutionDirectory.toString());
            //} else {
            //    edgeHostingAmounts = makeSolveLP(data, crf, numYears, captureTarget, edgeHostingAmounts, "");
            //}
        }
    }

    private static HashMap<Edge, Double> makeSolveLP(DataStorer data, double crf, double numYears, double captureTarget, HashMap<Edge, Double> edgeHostingAmounts, String solutionDirectory) {
        Source[] sources = data.getSources();
        Sink[] sinks = data.getSinks();
        LinearComponent[] linearComponents = data.getLinearComponents();
        int[] graphVertices = data.getGraphVertices();
        HashMap<Integer, HashSet<Integer>> neighbors = data.getGraphNeighbors();
        HashMap<Edge, Double> edgeConstructionCosts = data.getGraphEdgeConstructionCosts();

        //HashMap<Edge, Double> edgeRightOfWayCosts = data.getGraphEdgeRightOfWayCosts();
        HashMap<Source, Integer> sourceCellToIndex = new HashMap<>();
        HashMap<Integer, Source> sourceIndexToCell = new HashMap<>();
        HashMap<Sink, Integer> sinkCellToIndex = new HashMap<>();
        HashMap<Integer, Sink> sinkIndexToCell = new HashMap<>();
        HashMap<Integer, Integer> vertexCellToIndex = new HashMap<>();
        HashMap<Integer, Integer> vertexIndexToCell = new HashMap<>();
        HashMap<UnidirEdge, Integer> edgeToIndex = new HashMap<>();
        HashMap<Integer, UnidirEdge> edgeIndexToEdge = new HashMap<>();
        HashSet<Integer> sourceCells = new HashSet<>();
        HashSet<Integer> sinkCells = new HashSet<>();

        // Initialize cell/index maps
        for (int i = 0; i < sources.length; i++) {
            sourceCellToIndex.put(sources[i], i);
            sourceIndexToCell.put(i, sources[i]);
            sourceCells.add(sources[i].getCellNum());
        }
        for (int i = 0; i < sinks.length; i++) {
            sinkCellToIndex.put(sinks[i], i);
            sinkIndexToCell.put(i, sinks[i]);
            sinkCells.add(sinks[i].getCellNum());
        }
        for (int i = 0; i < graphVertices.length; i++) {
            vertexCellToIndex.put(graphVertices[i], i);
            vertexIndexToCell.put(i, graphVertices[i]);
        }
        int index = 0;
        for (Edge e : edgeConstructionCosts.keySet()) {
            UnidirEdge e1 = new UnidirEdge(e.v1, e.v2);
            edgeToIndex.put(e1, index);
            edgeIndexToEdge.put(index, e1);
            index++;

            UnidirEdge e2 = new UnidirEdge(e.v2, e.v1);
            edgeToIndex.put(e2, index);
            edgeIndexToEdge.put(index, e2);
            index++;
        }

        try {
            IloCplex cplex = new IloCplex();

            // variable: a
            IloNumVar[] a = new IloNumVar[sources.length];
            for (int s = 0; s < sources.length; s++) {
                a[s] = cplex.numVar(0, Double.MAX_VALUE, "a[" + s + "]");
            }

            // variable: b
            IloNumVar[] b = new IloNumVar[sinks.length];
            for (int s = 0; s < sinks.length; s++) {
                b[s] = cplex.numVar(0, Double.MAX_VALUE, "b[" + s + "]");
            }

            // variable: p
            IloNumVar[][] p = new IloNumVar[edgeToIndex.size()][linearComponents.length];
            for (int e = 0; e < edgeToIndex.size(); e++) {
                for (int c = 0; c < linearComponents.length; c++) {
                    p[e][c] = cplex.numVar(0, Double.MAX_VALUE, "p[" + e + "][" + c + "]");
                }
            }

            // constraint A: pipe capacity
            for (int e = 0; e < edgeToIndex.size(); e++) {
                for (int c = 0; c < linearComponents.length; c++) {
                    IloLinearNumExpr expr = cplex.linearNumExpr();
                    expr.addTerm(p[e][c], 1.0);
                    cplex.addLe(expr, linearComponents[c].getMaxCapacity());
                }
            }

            // constraint B: conservation of flow
            for (int src : graphVertices) {
                IloLinearNumExpr expr = cplex.linearNumExpr();
                for (int dest : neighbors.get(src)) {
                    UnidirEdge edge = new UnidirEdge(src, dest);
                    for (int c = 0; c < linearComponents.length; c++) {
                        expr.addTerm(p[edgeToIndex.get(edge)][c], 1.0);
                    }
                }

                for (int dest : neighbors.get(src)) {
                    UnidirEdge edge = new UnidirEdge(dest, src);
                    for (int c = 0; c < linearComponents.length; c++) {
                        expr.addTerm(p[edgeToIndex.get(edge)][c], -1.0);
                    }
                }

                // Set right hand side
                if (sourceCells.contains(src)) {
                    for (Source source : sources) {
                        if (source.getCellNum() == src) {
                            expr.addTerm(a[sourceCellToIndex.get(source)], -1.0);
                        }
                    }
                }
                if (sinkCells.contains(src)) {
                    for (Sink sink : sinks) {
                        if (sink.getCellNum() == src) {
                            expr.addTerm(b[sinkCellToIndex.get(sink)], 1.0);
                        }
                    }
                }
                cplex.addEq(expr, 0.0);
            }

            // constraint C: capture limit
            for (int s = 0; s < sources.length; s++) {
                IloLinearNumExpr expr = cplex.linearNumExpr();
                expr.addTerm(a[s], 1.0);
                cplex.addLe(expr, sources[s].getProductionRate());
            }

            // constraint D: injection limit
            for (int s = 0; s < sinks.length; s++) {
                IloLinearNumExpr expr = cplex.linearNumExpr();
                expr.addTerm(b[s], 1.0);
                cplex.addLe(expr, sinks[s].getCapacity() / numYears);
            }

            // constraint E: capture target
            IloLinearNumExpr expr = cplex.linearNumExpr();
            for (int s = 0; s < sources.length; s++) {
                expr.addTerm(a[s], 1.0);
            }
            cplex.addGe(expr, captureTarget);
            
            // constraint H: hardcoded values
            IloNumVar captureTargetVar = cplex.numVar(captureTarget, captureTarget, "captureTarget");
            IloLinearNumExpr h1 = cplex.linearNumExpr();
            h1.addTerm(captureTargetVar, 1.0);
            cplex.addEq(h1, captureTarget);
            
            IloNumVar crfVar = cplex.numVar(crf, crf, "crf");
            IloLinearNumExpr h2 = cplex.linearNumExpr();
            h2.addTerm(crfVar, 1.0);
            cplex.addEq(h2, crf);
            
            IloNumVar projectLengthVar = cplex.numVar(numYears, numYears, "projectLength");
            IloLinearNumExpr h3 = cplex.linearNumExpr();
            h3.addTerm(projectLengthVar, 1.0);
            cplex.addEq(h3, numYears);

            // objective
            IloLinearNumExpr objExpr = cplex.linearNumExpr();
            for (int s = 0; s < sources.length; s++) {
                objExpr.addTerm(a[s], sources[s].getCaptureCost());
            }

            for (int s = 0; s < sinks.length; s++) {
                objExpr.addTerm(b[s], sinks[s].getInjectionCost());
            }

            for (int e = 0; e < edgeToIndex.size(); e++) {
                for (int c = 0; c < linearComponents.length; c++) {
                    UnidirEdge unidirEdge = edgeIndexToEdge.get(e);
                    Edge bidirEdge = new Edge(unidirEdge.v1, unidirEdge.v2);

                    double fixedCost = (linearComponents[c].getConIntercept() * edgeConstructionCosts.get(bidirEdge)) * crf;
                    double variableCost = (linearComponents[c].getConSlope() * edgeConstructionCosts.get(bidirEdge)) * crf;
                    double coefficient = variableCost + (fixedCost / 1);
                    if (edgeHostingAmounts.containsKey(bidirEdge)) {
                        coefficient = variableCost + (fixedCost / edgeHostingAmounts.get(bidirEdge));
                    }
                    objExpr.addTerm(p[e][c], coefficient);
                }
            }

            // objective:
            IloObjective obj = cplex.minimize(objExpr);
            cplex.add(obj);
            cplex.setOut(null);

            // solve
            if (cplex.solve()) {
                for (int e = 0; e < edgeToIndex.size(); e++) {
                    for (int c = 0; c < linearComponents.length; c++) {
                        if (cplex.getValue(p[e][c]) > 0.0001) {
                            UnidirEdge unidirEdge = edgeIndexToEdge.get(e);
                            Edge bidirEdge = new Edge(unidirEdge.v1, unidirEdge.v2);
                            edgeHostingAmounts.put(bidirEdge, cplex.getValue(p[e][c]));
                        }
                    }
                }
            } else {
                System.out.println("Not Feasible");
            }

            if (solutionDirectory != null && solutionDirectory != "") {
                cplex.exportModel(solutionDirectory + "/flowCap.mps");
                cplex.writeSolution(solutionDirectory + "/solution.sol");
                Solution soln = DataInOut.loadSolution(solutionDirectory, -1);
                System.out.println("Solution Cost = " + soln.getTotalCost());
            }
            cplex.clearModel();
        } catch (IloException ex) {
            Logger.getLogger(FlowHeuristic.class.getName()).log(Level.SEVERE, null, ex);
        }

        return edgeHostingAmounts;
    }
}
