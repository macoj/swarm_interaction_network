package pso;
/**
*
* @author Bruno Andrade
*/

import java.io.File;
import java.io.FileNotFoundException;
import java.io.PrintWriter;
import java.util.logging.Level;
import java.util.logging.Logger;

public class Trigger {

	static String[] arg = new String[1];

	public static void main(String[] args){
		String pathSeparator = File.separator;
		File rootPath = new File(System.getProperty("user.dir") + pathSeparator + "psoResults");
        if (!rootPath.exists()) {
            rootPath.mkdir();
        }
		
		for (String function : Holder.functions) {
			File functionFolder = new File(rootPath + pathSeparator + "F" + function);
			if (!functionFolder.exists()) {
				functionFolder.mkdir();
			}
			for (PSO.TOPOLOGY topology : Holder.topologies) {
				File topologyFolder = new File(rootPath + pathSeparator + "F" + function + pathSeparator + topology.toString());
				if (!topologyFolder.exists()) {
					topologyFolder.mkdir();
				}
				for (PSO.TOPOLOGY_MECHANISM topology_mechanism : Holder.topologyMechanims) {
					for(int i=0; i<Holder.EXECUTIONS; i++) {
						try {
							arg[0] = function;
							PSO pso = new PSO(arg);
							pso.setTopology(topology);
							pso.setTopologyMechanism(topology_mechanism);
							String filename = "psoResults_F" + function +"_"+topology;
							if (pso.swarm_topology_mechanism != null) {
								filename += "_"+ pso.swarm_topology_mechanism.toString();
							}   
							filename += "_" + (i+1) +".txt";
							File psoFile = new File(rootPath + pathSeparator + "F" + function + pathSeparator + topology.toString() + pathSeparator +  pathSeparator + filename);
							System.out.println("File being created at " + psoFile.getAbsolutePath());
							PrintWriter printWriter = new PrintWriter(psoFile);
							pso.setPrinter(printWriter);
							Thread thread = new Thread(pso);
							thread.start();
						} catch (FileNotFoundException ex) {
							Logger.getLogger(Trigger.class.getName()).log(Level.SEVERE, null, ex);
						} finally {

						}
					}
				}
			}
		}   
	}
}

