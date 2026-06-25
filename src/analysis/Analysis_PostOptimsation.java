package analysis;

import java.io.File;
import java.io.FileFilter;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.Map.Entry;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import optimisation.Abstract_Optimisation;
import util.StaticMethods;

public class Analysis_PostOptimsation {
	
	public final static String filenameformat_optResultsByMap = "OptResult_%d.csv";
	
	private final static String filenameformat_bestResByMap = "OptResult_BestPerMap.csv";	
	private final static Pattern pattern_param_list = Pattern.compile(
			String.format(Abstract_Optimisation.fileformat_opt_outcomes,"(.*)_(\\d+)\\")						
//			"OptProgress_ParamList_(.*)_(\\d+)\\.csv"			
			);
	
	private static class OptResult_Ent {
		final Double residue;
		final String line;
		final String fName;
		final Long cMap;
		
		private OptResult_Ent(String line, String fName) {
			this.line = line;
			this.fName = fName;
			String[] lineEnt = line.split(",");
			cMap = Long.parseLong(lineEnt[0]);
			residue = Double.valueOf(lineEnt[lineEnt.length - 1]);
		}
	}

	public static void analysis(String baseSimDirName, String csvDirName) throws FileNotFoundException, IOException {

		File baseDir = new File(baseSimDirName);

		File[] optBestDirs = baseDir.listFiles(new FileFilter() {
			@Override
			public boolean accept(File pathname) {
				return pathname.isDirectory() && new File(pathname, csvDirName).isDirectory();
			}
		});

		Arrays.sort(optBestDirs);

		HashMap<String, ArrayList<StringBuilder>> map_lines_collections = new HashMap<>();

		for (File optBestDir : optBestDirs) {

			File csv_dir = new File(optBestDir, csvDirName);

			File[] csv_list = csv_dir.listFiles(new FileFilter() {
				@Override
				public boolean accept(File pathname) {
					return pathname.getName().endsWith("csv");
				}
			});
			for (File csv : csv_list) {
				System.out.printf("Looking into %s\n", csv.getAbsolutePath());
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
		File outputDir = new File(baseDir, csvDirName);
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

		Analysis_PostOptimsation.extractedOptResultsByMap(new File[] { baseDir });

	}

	public static void extractedOptResultsByMap(File[] resultDirs_list) throws FileNotFoundException, IOException {

		HashMap<Long, ArrayList<OptResult_Ent>> optResultByCMap = new HashMap<>();

		Comparator<OptResult_Ent> cmp = new Comparator<OptResult_Ent>() {
			@Override
			public int compare(OptResult_Ent o1, OptResult_Ent o2) {
				return o1.residue.compareTo(o2.residue);
			}
		};

		for (File resultDirs : resultDirs_list) {
			File[] result_files = resultDirs.listFiles(new FileFilter() {
				@Override
				public boolean accept(File pathname) {
					return pattern_param_list.matcher(pathname.getName()).matches();
				}
			});

			File[] baseSeedList = resultDirs.listFiles(new FileFilter() {
				@Override
				public boolean accept(File pathname) {
					return pathname.isDirectory() && pathname.getName().startsWith("Seed_List_");
				}
			});

			ArrayList<String> noResult_baseSeedDirNames = new ArrayList<>(baseSeedList.length);
			for (File f : baseSeedList) {
				noResult_baseSeedDirNames.add(f.getName());
			}
			Collections.sort(noResult_baseSeedDirNames);

			String header = null;

			HashMap<Integer, ArrayList<String>> num_result_by_seed_dir = new HashMap<>();

			for (File f : result_files) {
				Matcher m = pattern_param_list.matcher(f.getName());
				m.matches();
				String seed_dir = m.group(1);

				int pt = Collections.binarySearch(noResult_baseSeedDirNames, seed_dir);
				if (pt >= 0) {
					noResult_baseSeedDirNames.remove(pt);
				}

				String[] lines = util.Util_7Z_CSV_Entry_Extract_Callable.extracted_lines_from_text(f);
				if (header == null) {
					header = lines[0];
				}

				ArrayList<String> numResultsArr = num_result_by_seed_dir.get(lines.length - 1);
				if (numResultsArr == null) {
					numResultsArr = new ArrayList<>();
					num_result_by_seed_dir.put(lines.length - 1, numResultsArr);
				}
				numResultsArr.add(seed_dir);

				for (int i = 1; i < lines.length; i++) {
					OptResult_Ent arrayEnt = new OptResult_Ent(lines[i], f.getName());
					Long cMap_seed = arrayEnt.cMap;

					ArrayList<OptResult_Ent> optResultArr = optResultByCMap.get(cMap_seed);
					if (optResultArr == null) {
						optResultArr = new ArrayList<>();
						optResultByCMap.put(cMap_seed, optResultArr);
					}

					int insertPt = Collections.binarySearch(optResultArr, arrayEnt, cmp);
					if (insertPt < 0) {
						insertPt = ~insertPt;
					}
					optResultArr.add(insertPt, arrayEnt);
				}
				System.out.printf("Reading of %s completed. # entries = %d\n", f.getName(), lines.length - 1);
			}

			PrintWriter pri_best_per_map = new PrintWriter(new File(resultDirs, filenameformat_bestResByMap));
			pri_best_per_map.printf("DIR_NAME,%s\n", header);

			for (Entry<Long, ArrayList<OptResult_Ent>> ent : optResultByCMap.entrySet()) {
				File optResultCMap_CSV = new File(resultDirs, String.format(filenameformat_optResultsByMap, ent.getKey()));
				PrintWriter pWri = new PrintWriter(optResultCMap_CSV);
				pWri.printf("DIR_NAME,%s\n", header);

				ArrayList<Long> bestRowBySimSeed = new ArrayList<>();

				for (OptResult_Ent arrEnt : ent.getValue()) {

					Long sim_seed = Long.valueOf(arrEnt.line.split(",")[1]);

					int pt = Collections.binarySearch(bestRowBySimSeed, sim_seed);

					if (pt < 0) {
						pri_best_per_map.printf("%s,%s\n", arrEnt.fName, arrEnt.line);
						bestRowBySimSeed.add(~pt, sim_seed);
					}

					pWri.printf("%s,%s\n", arrEnt.fName, arrEnt.line);
				}
				pWri.close();
				System.out.printf("Result of CMAP=%d generated at %s\n", ent.getKey(), optResultCMap_CSV.getName());
			}
			pri_best_per_map.close();

			// Prioritise next seed list based on those haven't run yet

			PrintWriter pri_dir_list = new PrintWriter(new File(resultDirs, "SeedDirListing_Prior.txt"));

			StringBuilder info = new StringBuilder();

			info.append(String.format("# seed list directories in total = %d\n", baseSeedList.length));
			info.append(String.format("# directories with result = %d\n",
					baseSeedList.length - noResult_baseSeedDirNames.size()));
			info.append(String.format("# directories with no result yet = %d\n", noResult_baseSeedDirNames.size()));

			Collections.shuffle(noResult_baseSeedDirNames);
			for (String dirName : noResult_baseSeedDirNames) {
				pri_dir_list.printf("   %s\n", dirName);
			}

			// Generate list for next set of result
			Integer[] res_arr = num_result_by_seed_dir.keySet().toArray(new Integer[0]);
			Arrays.sort(res_arr);
			for (Integer nRes : res_arr) {
				ArrayList<String> dirNames = num_result_by_seed_dir.get(nRes);
				Collections.shuffle(dirNames);
				info.append(String.format("# directories with %d results = %d\n", nRes, dirNames.size()));
				for (String dirName : dirNames) {
					pri_dir_list.printf("   %s\n", dirName);
				}
			}

			pri_dir_list.print("\n");
			pri_dir_list.print(info.toString());
			pri_dir_list.close();

			System.out.println(info.toString());

		}
	}

}
