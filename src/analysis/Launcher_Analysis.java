package analysis;

import java.io.File;
import java.io.FileFilter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Map.Entry;

import optimisation.MenB_RMP_NM_Optimistion;
import sim.Simulation_ClusterModelTransmission;
import util.StaticMethods;

public class Launcher_Analysis {

	public static void main(String[] args) throws IOException {

		String usageInfo = "Usage: java -jar Analyse_MenB_Vaccine.jar BASEDIR_SIM PATH_REGION_MAPPING PATH_GRP_SIZE <-flag=BINARY_FLAG_OPTIONS> <-printProgress=TF>"
				+ "\n  or java -jar Analyse_MenB_Vaccine.jar -opt BASEDIR_SIM SEED_DIR_NAME <-printProgress=TF> <-optType=OPT_TYPE>"
				+ "\n  or java -jar Analyse_MenB_Vaccine.jar -optBestFit BASEDIR_SIM CSV_FOLDER_NAME";

		boolean correctUsage = false;
		if (args.length > 0) {

			if ("-opt".equals(args[0]) && args.length >= 3) {
				// Optimisation
				correctUsage = true;
				MenB_RMP_NM_Optimistion opt = new MenB_RMP_NM_Optimistion(args[1], args[2]);
				for (int i = 2; i < args.length; i++) {
					if (args[i].startsWith(Simulation_ClusterModelTransmission.LAUNCH_ARGS_PRINT_PROGRESS)) {
						opt.setObjFunc_PrintProgress(Boolean.parseBoolean(args[i].split("=")[1]));
					}
					if (args[i].startsWith("-optType")) {
						opt.setOptType(Integer.parseInt(args[i].split("=")[1]));
					}
				}
				opt.runOptimisation();
			} else if ("-optBestFit".equals(args[0]) && args.length >= 3) {
				correctUsage = true;
				File baseDir = new File(args[1]);
				File[] optBestDirs = baseDir.listFiles(new FileFilter() {
					@Override
					public boolean accept(File pathname) {
						return pathname.isDirectory() && new File(pathname, args[2]).isDirectory();
					}
				});

				Arrays.sort(optBestDirs);

				HashMap<String, ArrayList<StringBuilder>> map_lines_collections = new HashMap<>();

				for (File optBestDir : optBestDirs) {

					File csv_dir = new File(optBestDir, args[2]);

					File[] csv_list = csv_dir.listFiles(new FileFilter() {
						@Override
						public boolean accept(File pathname) {
							return pathname.getName().endsWith("csv");
						}
					});
					for (File csv : csv_list) {
						String key = csv.getName();
						String[] lines = StaticMethods.extracted_lines_from_text(csv);
						StringBuilder bld;
						String[] ent;

						ArrayList<StringBuilder> lines_arr = map_lines_collections.get(key);

						if (lines_arr == null) {
							lines_arr = new ArrayList<>();
							map_lines_collections.put(key, lines_arr);
						}
						if (lines_arr.size() == 0) {
							// Add header and first time column
							for (String s : lines) {
								bld = new StringBuilder();
								bld.append(s.split(",")[0]);
								lines_arr.add(bld);
							}
						}

						// Header
						bld = lines_arr.get(0);
						ent = lines[0].split(",");

						for (int c = 1; c < ent.length; c++) {
							bld.append(',');
							bld.append(String.format("%s:%s", optBestDir.getName(), ent[c]));
						}

						// Entry
						for (int r = 1; r < lines_arr.size(); r++) {
							bld = lines_arr.get(r);
							if (r < lines.length) {
								ent = lines[r].split(",");
								for (int c = 1; c < ent.length; c++) {
									bld.append(',');
									bld.append(ent[c]);
								}
							} else {
								System.err.printf("Error! Number of lines in %s doesn't match with first result.\n");
								System.exit(-1);
							}
						}
					}
				}

				// Print grouped result
				File outputDir = new File(baseDir, args[2]);
				outputDir.mkdirs();

				for (Entry<String, ArrayList<StringBuilder>> ent : map_lines_collections.entrySet()) {

					File outputCSV = new File(outputDir, ent.getKey());
					PrintWriter pWri = new PrintWriter(outputCSV);
					ArrayList<StringBuilder> lines_arr = ent.getValue();

					for (StringBuilder bld : lines_arr) {
						pWri.println(bld.toString());
					}

					pWri.close();
					System.out.printf("Consolidate CSV %s generated.\n", outputCSV.getAbsolutePath());

				}

			} else if (args.length >= 3) {
				// Time trend
				correctUsage = true;
				Analysis_PostSim_ExtractTimeTrends analyse = new Analysis_PostSim_ExtractTimeTrends();
				analyse.analyse(args);
			}
		}

		if (!correctUsage) {
			System.out.println(usageInfo);
			System.exit(0);
		}

	}

}
