package analysis;

import java.io.File;
import java.io.FileFilter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.HashMap;
import java.util.Map;
import java.util.Map.Entry;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import sim.Simulation_ClusterModelTransmission;
import util.StaticMethods;

public class Analysis_PostSim_ExtractTimeTrends {

	Pattern pattern_sim = Pattern.compile("Seed_List.*_(\\d+)");

	// Check completeness
	Pattern[] pattern_check_file = new Pattern[] { //
			// Pattern.compile("Incidence_Person_.*csv\\.7z"), //
			// Pattern.compile("Infectious_by_GrpLoc_.*csv\\.7z"), //
			// Pattern.compile("Vaccine_Coverage_.*csv\\.7z"), //
			// Pattern.compile("Treatment_Person_-?\\d+.*csv\\.7z"), //
			// Pattern.compile("InfectHist_.*csv\\.7z")
			Pattern.compile("Treatment_by_GrpLoc_.*csv\\.7z"), Pattern.compile("Infectious_by_GrpRegion_.*csv\\.7z"), };

	// Extract time trend
	Map<String, Pattern> map_timetrend = Map.ofEntries( //
			Map.entry("Timetrend_Incidence_Person_%s.csv", Pattern.compile("Incidence_Person_.*csv\\.7z")), //
			Map.entry("Timetrend_Treatment_Person_%s.csv", Pattern.compile("Treatment_Person_.*csv\\.7z")), //
			Map.entry("Timetrend_Treatment_by_GrpLoc_%s.csv", Pattern.compile("Treatment_by_GrpLoc_.*csv\\.7z")), //
			Map.entry("Timetrend_Infectious_by_GrpLoc_%s.csv", Pattern.compile("Infectious_by_GrpLoc_.*csv\\.7z")), //
			Map.entry("Timetrend_Infectious_by_GrpRegion_%s.csv", Pattern.compile("Infectious_by_GrpRegion_.*csv\\.7z")) //
	);

	String[][] grp_extract_array = new String[][] { //
			new String[] { //
					"Timetrend_Treatment_by_GrpLoc_%s.csv", //
					"Timetrend_Treatment_by_GrpLoc_%s.csv", //
			}, //
			new String[] { //
					"Timetrend_Treatment_by_GrpLoc_Grp_%d_Loc_%d.csv", //
					"Timetrend_Treatment_by_GrpLoc_Grp_%d_Loc_%d.csv", //
			}, //
			new String[] { //
					"[56789]", // Indigenous female
					"1[56789]", // Non-Indigenous female
			}//

	};

	String[][] region_extract_array = new String[][] { //
			new String[] { //
					"Timetrend_Treatment_by_GrpLoc_%s.csv", //
					"Timetrend_Infectious_by_GrpLoc_%s.csv", //

			}, //
			new String[] { //
					"Timetrend_Treatment_by_GrpLoc_Grp_%d_Loc_%d.csv", //
					"Timetrend_Infectious_by_GrpLoc_Grp_%d_Loc_%d.csv", //
			}, //
			new String[] { //
					"Timetrend_Treatment_%s_at_%s.csv", //
					"Timetrend_Infectious_%s_at_%s.csv", //
			} //
	};

	public void setMap_timetrend(Map<String, Pattern> map_timetrend) {
		this.map_timetrend = map_timetrend;
	}

	public void analyse(String[] args) throws IOException {

		File basedir_sim = new File(args[0]);
		File file_regions_mapping = new File(args[1]);
		File file_grp_size = new File(args[2]);

		boolean flag_check_completeness = true; // 1
		boolean flag_extract_timetrend = true; // 2
		boolean flag_grp_region_extract = true; // 4
		boolean flag_infection_hist_pid = true; // 8

		boolean flag_printProgress = false;
		boolean keepGrpRegionCSVs = false;

		for (int i = 3; i < args.length; i++) {
			if (args[i].startsWith("-flag=")) {
				int flag = Integer.parseInt(args[i].substring("-flag=".length()));
				flag_check_completeness = ((1 << 0 & flag) != 0);
				flag_extract_timetrend = ((1 << 1 & flag) != 0);
				flag_grp_region_extract = ((1 << 2 & flag) != 0);
				flag_infection_hist_pid = ((1 << 3 & flag) != 0);
			}
			if (args[i].startsWith("-keepGrpRegionCSV=")) {
				keepGrpRegionCSVs = Boolean.valueOf(args[i].substring("-flag=".length()));
			}

			if (args[i].startsWith(Simulation_ClusterModelTransmission.LAUNCH_ARGS_PRINT_PROGRESS)) {
				String[] ent = args[i].split("=");
				flag_printProgress = Boolean.parseBoolean(ent[ent.length - 1]);
			}
		}

		Map<String, File> map_timetrend_dir = Map.ofEntries( //
				Map.entry("Timetrend_Incidence_Person_%s.csv", new File(basedir_sim, "Timetrend_Incidence")) //
				, Map.entry("Timetrend_Treatment_Person_%s.csv", new File(basedir_sim, "Timetrend_Treatment_Person")) //
				, Map.entry("Timetrend_Treatment_by_GrpLoc_%s.csv", new File(basedir_sim, "Timetrend_Treatment_GrpLoc")) //
				,
				Map.entry("Timetrend_Infectious_by_GrpLoc_%s.csv", new File(basedir_sim, "Timetrend_Infectious_GrpLoc")) //
				, Map.entry("Timetrend_Infectious_by_GrpRegion_%s.csv",
						new File(basedir_sim, "Timetrend_Infectious_GrpReg")) //
		);

		int model_pop_size = 1000000;

		String[] datasel_by_grp_key = grp_extract_array[0];
		String[] datasel_by_grp_csv_strformat = grp_extract_array[1];
		String[] datasel_by_grp_output_regEx = grp_extract_array[2];

		// Regional mapping
		String[] datasel_by_regions_key = region_extract_array[0];
		String[] datasel_by_regions_csv_strformat = region_extract_array[1];
		String[] datasel_by_regions_output_fname = region_extract_array[2];

		Map<String, int[]> map_timetrend_grpIncl = Map.ofEntries( //
				Map.entry("Indigenous", new int[] { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9 }) // Indigenous
				, Map.entry("Non-Indigenous", new int[] { 10, 11, 12, 13, 14, 15, 16, 17, 18, 19 }) // Non-Indigenous
		);

		// Code start
		if (flag_printProgress) {
			System.out.printf("Checking/Analysing %s:\n", basedir_sim.getAbsolutePath());
		}

		File[] dirs_sim = basedir_sim.listFiles(new FileFilter() {
			@Override
			public boolean accept(File pathname) {
				return pathname.isDirectory() && pattern_sim.matcher(pathname.getName()).matches();
			}
		});

		// Check completeness of results
		if (flag_check_completeness) {
			ArrayList<File> incomplete_dirs = new ArrayList<>();
			for (File simDir : dirs_sim) {
				File[] checkFile = simDir.listFiles(new FileFilter() {
					@Override
					public boolean accept(File pathname) {
						boolean res = false;
						for (Pattern p : pattern_check_file) {
							res |= p.matcher(pathname.getName()).matches();
						}
						return res;
					}
				});

				if (checkFile.length != pattern_check_file.length) {
					incomplete_dirs.add(simDir);
				}

			}

			if (flag_printProgress) {
				System.out.printf("# of incomplete sim directory = %d\n", incomplete_dirs.size());
			}

			if (!incomplete_dirs.isEmpty()) {

				PrintWriter pri_batch_del = new PrintWriter(new File(basedir_sim, "batch_del")) {
					@Override
					public void println() {
						write('\n');
					}
				};

				pri_batch_del.println("#!/bin/bash");
				pri_batch_del.println("echo \"Batch deleting...\"");

				System.out.printf("Incomplete directories (%d in total):\n", incomplete_dirs.size());
				for (File incompleteDir : incomplete_dirs) {
					System.out.printf("   -seedMap=%s/%s.csv\n", incompleteDir.getName(), incompleteDir.getName());
					pri_batch_del.printf("find ./%s -type f \\! -name \"Seed_List*csv\" -delete\n",
							incompleteDir.getName());

				}

				pri_batch_del.println("echo \"Batch deleting completed\"");

				pri_batch_del.close();				
				System.exit(0);
			}
		}

		if (flag_extract_timetrend) {
			Pattern pattern_zip_ent = Pattern.compile("\\[(.*),(\\d+)\\].*");

			for (Entry<String, Pattern> timetrendEntry : map_timetrend.entrySet()) {
				Pattern tt_pattern = timetrendEntry.getValue();
				String tt_name = timetrendEntry.getKey();

				// K = ColName, V = Map< DirName:SeedFileName:SeedFileRow, Entries>
				HashMap<String, HashMap<String, String[]>> timetrendMap = new HashMap<>();

				ArrayList<String> error_dirs = new ArrayList<>();

				for (File simDir : dirs_sim) {
					File[] sel_file = simDir.listFiles(new FileFilter() {
						@Override
						public boolean accept(File pathname) {
							return tt_pattern.matcher(pathname.getName()).matches();
						}
					});
					if (sel_file.length == 1) {
						if (flag_printProgress) {
							System.out.printf("Extracting data from %s\n", sel_file[0].getAbsolutePath());
						}
						try {
							HashMap<String, ArrayList<String[]>> lineMap;

							if (sel_file[0].getName().endsWith("7z")) {
								lineMap = util.Util_7Z_CSV_Entry_Extract_Callable.extractedLinesFrom7Zip(sel_file[0]);
							} else {
								// Direct CSV version
								lineMap = new HashMap<>();
								String[] lines = util.Util_7Z_CSV_Entry_Extract_Callable
										.extracted_lines_from_text(sel_file[0]);
								ArrayList<String[]> lineSp = new ArrayList<>();
								for (String line : lines) {
									lineSp.add(line.split(","));
								}
								lineMap.put(sel_file[0].getName(), lineSp);
							}

							for (Entry<String, ArrayList<String[]>> lineEnt : lineMap.entrySet()) {
								Matcher m_filename = pattern_zip_ent.matcher(lineEnt.getKey());
								if (m_filename.matches()) {
									String map_key = String.format("%s:%s:%s", simDir.getName(), m_filename.group(1),
											m_filename.group(2));

									String[] col_name = lineEnt.getValue().remove(0);
									for (int c = 0; c < col_name.length; c++) {
										if (!col_name[c].endsWith("_null")) { // Skip null location
											HashMap<String, String[]> timetrendMapByCol = timetrendMap.get(col_name[c]);
											if (timetrendMapByCol == null) {
												timetrendMapByCol = new HashMap<>();
												timetrendMap.put(col_name[c], timetrendMapByCol);
											}
											if (lineEnt.getValue().size() > 0) {
												String[] ent = new String[lineEnt.getValue().size()];
												for (int i = 0; i < lineEnt.getValue().size(); i++) {
													ent[i] = lineEnt.getValue().get(i)[c];
												}
												timetrendMapByCol.put(map_key, ent);
											}
										}
									}

								} else {
									System.err.printf("Error! Unexpected file name %s in %s\n", lineEnt.getKey(),
											sel_file[0].getAbsolutePath());
								}

							}
						} catch (Exception ex) {
							System.err.printf("Error! Following %s encountered during reading of %s entry.\n",
									ex.getMessage(), sel_file[0].getAbsolutePath());
							ex.printStackTrace(System.err);
							error_dirs.add(simDir.getName());

						}

					}

				}

				if (!error_dirs.isEmpty()) {
					System.out.printf("The following %s sim directory has error(s):\n", error_dirs.size());
					for (String err_dir : error_dirs) {
						System.out.printf("   -seedMap=%s/%s.csv\n", err_dir, err_dir);
					}
				}

				if (timetrendMap.isEmpty()) {
					System.out.printf("Extract of time trend for \"%s\" skipped.\n", tt_name);
				} else {
					// Printing of outputs
					String[] time_col = new String[0];
					ArrayList<String> sim_key_array = new ArrayList<>();

					HashMap<String, String[]> map_time = timetrendMap.remove("Time");

					for (Entry<String, String[]> ent : map_time.entrySet()) {
						sim_key_array.add(ent.getKey());
						String[] tCol = ent.getValue();
						if (tCol.length > time_col.length) {
							time_col = tCol;
						}
					}

					sim_key_array.sort(new Comparator<String>() {
						@Override
						public int compare(String o1, String o2) {
							String[] s1 = o1.split(":");
							String[] s2 = o2.split(":");
							int res = 0;

							for (int p = 0; p < Math.min(s1.length, s2.length) && res == 0; p++) {
								switch (p) {
								case 2:
									res = Integer.compare(Integer.parseInt(s1[p]), Integer.parseInt(s2[p]));
									break;
								default:
									res = s1[p].compareTo(s2[p]);
								}
							}

							return res;
						}
					});

					for (Entry<String, HashMap<String, String[]>> tt_map_by_col : timetrendMap.entrySet()) {
						String col_name = tt_map_by_col.getKey();
						HashMap<String, String[]> tt_map_by_sim = tt_map_by_col.getValue();

						StringBuilder[] linesbuilder = new StringBuilder[time_col.length];
						for (int r = 0; r < linesbuilder.length; r++) {
							linesbuilder[r] = new StringBuilder();
							linesbuilder[r].append(time_col[r]);
						}

						for (String sim_key : sim_key_array) {
							String[] ent = tt_map_by_sim.get(sim_key);
							for (int r = 0; r < linesbuilder.length; r++) {
								linesbuilder[r].append(',');
								if (r < ent.length) {
									linesbuilder[r].append(ent[r]);
								} else {
									linesbuilder[r].append("NaN");
								}

							}
						}

						File basedir = basedir_sim;

						if (map_timetrend_dir.containsKey(tt_name)) {
							basedir = map_timetrend_dir.get(tt_name);
							basedir.mkdirs();
						}

						File f_timetrend = new File(basedir, String.format(tt_name, col_name));
						PrintWriter p_f_timetrend = new PrintWriter(f_timetrend);

						p_f_timetrend.print(col_name);
						for (String sim_key : sim_key_array) {
							p_f_timetrend.print(',');
							p_f_timetrend.print(sim_key);
						}

						p_f_timetrend.println();

						for (int r = 0; r < linesbuilder.length; r++) {
							p_f_timetrend.println(linesbuilder[r].toString());
						}

						p_f_timetrend.close();

						if (flag_printProgress) {
							System.out.printf("Time trend for %s generated at %s\n", col_name,
									f_timetrend.getAbsolutePath());
						}

					}

					if (flag_printProgress) {
						System.out.printf("Time trend extract for %s completed.\n", timetrendEntry.getKey());
					}
				}
			}

		}

		if (flag_grp_region_extract) {

			for (int resPt = 0; resPt < datasel_by_grp_csv_strformat.length; resPt++) {
				File data_dir = basedir_sim;
				if (map_timetrend_dir.containsKey(datasel_by_grp_key[resPt])) {
					data_dir = map_timetrend_dir.get(datasel_by_grp_key[resPt]);
				}
				String csv_format_by_grp = datasel_by_grp_csv_strformat[resPt]
						.replaceFirst("%d", datasel_by_grp_output_regEx[resPt]).replaceAll("%d", "\\\\d+");
				Pattern pattern_csv = Pattern.compile(csv_format_by_grp);
				String fname_output = datasel_by_grp_csv_strformat[resPt]
						.replaceFirst("%d", datasel_by_grp_output_regEx[resPt]).replaceAll("%d", "ALL");

				StringBuilder header = new StringBuilder();
				header.append(fname_output);

				File[] src_files = data_dir.listFiles(new FileFilter() {
					@Override
					public boolean accept(File pathname) {
						return pattern_csv.matcher(pathname.getName()).matches();
					}
				});

				ArrayList<double[]> val = new ArrayList<>();
				double[] valEnt;

				for (File src_file : src_files) {
					String[] src_lines = StaticMethods.extracted_lines_from_text(src_file);
					String[] src_header = src_lines[0].split(",");

					if (val.size() == 0) {
						for (int c = 1; c < src_header.length; c++) {
							header.append(',');
							header.append(src_header[c]);
						}
						for (int r = 1; r < src_lines.length; r++) {
							valEnt = new double[src_header.length];
							Arrays.fill(valEnt, Double.NaN);
							val.add(valEnt);
						}
					}

					for (int r = 1; r < src_lines.length; r++) {
						valEnt = val.get(r - 1);
						String[] src_ent = src_lines[r].split(",");
						if (Double.isNaN(valEnt[0])) {
							Arrays.fill(valEnt, 0);
							valEnt[0] = Integer.parseInt(src_ent[0]);
						}
						for (int c = 1; c < src_ent.length; c++) {
							valEnt[c] += Double.parseDouble(src_ent[c]);
						}

					}

				}

				// Print output
				String[] lines = new String[val.size() + 1];
				lines[0] = header.toString();
				for (int i = 1; i < lines.length; i++) {
					StringBuilder ent = new StringBuilder();
					valEnt = val.get(i - 1);
					if (!Double.isNaN(valEnt[0])) {
						ent.append((int) valEnt[0]);
						for (int c = 1; c < valEnt.length; c++) {
							ent.append(',');
							ent.append(valEnt[c]);
						}
					}

					lines[i] = ent.toString();
				}

				File output_file = new File(data_dir, fname_output);
				PrintWriter pWri_output = new PrintWriter(output_file);
				for (String line : lines) {
					pWri_output.println(line);
				}
				pWri_output.close();
			}

			for (int resPt = 0; resPt < datasel_by_regions_csv_strformat.length; resPt++) {

				File data_dir = basedir_sim;
				if (map_timetrend_dir.containsKey(datasel_by_regions_key[resPt])) {
					data_dir = map_timetrend_dir.get(datasel_by_regions_key[resPt]);
				}

				// K = region ID, V = popSize by grp
				HashMap<Integer, int[]> regions_popSize_raw = new HashMap<>();
				int raw_popSize = 0;
				// K = region ID, V = Grp->{city,region,remote}
				HashMap<Integer, HashMap<Integer, double[]>> regions_weight = new HashMap<>();

				String[] regions_size = util.Util_7Z_CSV_Entry_Extract_Callable
						.extracted_lines_from_text(file_grp_size);
				for (int r = 1; r < regions_size.length; r++) {
					String[] ent = regions_size[r].split(",");
					Integer region_id = Integer.valueOf(ent[0]);
					int[] grp_size = new int[ent.length - 2];
					for (int g = 2; g < ent.length; g++) {
						int gSize = Integer.parseInt(ent[g]);
						grp_size[g - 2] = gSize;
						raw_popSize += gSize;
					}
					regions_popSize_raw.put(region_id, grp_size);
				}

				String[] regions_list = util.Util_7Z_CSV_Entry_Extract_Callable
						.extracted_lines_from_text(file_regions_mapping);

				String[] regions_headers = regions_list[0].split(",");

				for (int i = 1; i < regions_list.length; i++) {
					String[] ent = regions_list[i].split(",");
					Integer reg_id = Integer.valueOf(ent[0]);
					Integer grp_id = Integer.valueOf(ent[1]);
					double[] val = new double[ent.length - 2];

					double pop_sum = 0;

					// Fill val
					for (int c = 2; c < ent.length; c++) {
						double entVal = Double.parseDouble(ent[c]);
						val[c - 2] = entVal;
						pop_sum += val[c - 2];
					}

					// Adjust weight
					for (int c = 0; c < val.length; c++) {
						val[c] = val[c] / pop_sum;
					}

					HashMap<Integer, double[]> grp_ent = regions_weight.get(reg_id);

					if (grp_ent == null) {
						grp_ent = new HashMap<>();
						regions_weight.put(reg_id, grp_ent);
					}

					grp_ent.put(grp_id, val);

				}

				Integer[] region_ids = regions_weight.keySet().toArray(new Integer[0]);
				Arrays.sort(region_ids);

				File file_adj_pop_size_by_region = new File(basedir_sim, "Adj_pop_size.csv");
				if (!file_adj_pop_size_by_region.exists()) {

					PrintWriter pWri_adj_pop_size = new PrintWriter(file_adj_pop_size_by_region);
					pWri_adj_pop_size.printf(
							"Region_Id,Grp," + "CITY_Adj_to_%1$d,REGIONAL_Adj_to_%1$d,REMOTE_Adj_to_%1$d\n",
							model_pop_size);

					for (Integer region_id : region_ids) {
						HashMap<Integer, double[]> weight_by_grp = regions_weight.get(region_id);

						Integer[] grp_ids = weight_by_grp.keySet().toArray(new Integer[0]);
						Arrays.sort(grp_ids);

						int[] popSize_raw = regions_popSize_raw.get(region_id);

						for (Integer grp_id : grp_ids) {
							double[] region_weight = weight_by_grp.get(grp_id);

							// Print output
							pWri_adj_pop_size.print(region_id);
							pWri_adj_pop_size.print(',');
							pWri_adj_pop_size.print(grp_id);

							double adjPopVal = 0;
							adjPopVal = Math
									.round(((double) popSize_raw[grp_id.intValue()]) / raw_popSize * model_pop_size);

							for (int r = 0; r < region_weight.length; r++) {
								pWri_adj_pop_size.print(',');
								pWri_adj_pop_size.print(adjPopVal * region_weight[r]);

							}

							pWri_adj_pop_size.println();
						}
					}
					pWri_adj_pop_size.close();
					System.out.printf("Adjusted popSize by region exported to %s\n",
							file_adj_pop_size_by_region.getAbsolutePath());
				}

				String[] common_sim_header = null;

				Pattern pattern_datasel_csv = Pattern
						.compile(datasel_by_regions_csv_strformat[resPt].replace("%d", "(\\d+)"));

				for (Entry<String, int[]> ent : map_timetrend_grpIncl.entrySet()) {
					int[] grps = ent.getValue();
					File[] datasel_files = data_dir.listFiles(new FileFilter() {
						@Override
						public boolean accept(File pathname) {
							boolean res;
							String fileName = pathname.getName();
							Matcher m = pattern_datasel_csv.matcher(fileName);
							res = m.matches();
							if (res) {
								int grpNum = Integer.parseInt(m.group(1));
								res = (Arrays.binarySearch(grps, grpNum)) >= 0;
							}
							return res;
						}
					});

					// double[ra_index][time][sim_num]
					double[][][] val_all = null;

					for (File datasel_file : datasel_files) {
						Matcher m = pattern_datasel_csv.matcher(datasel_file.getName());
						if (m.matches()) {
							String[] lines = util.Util_7Z_CSV_Entry_Extract_Callable
									.extracted_lines_from_text(datasel_file);

							String[] header = lines[0].split(",");
							Integer region_id = Integer.valueOf(m.group(2));
							Integer grp_id = Integer.valueOf(m.group(1));

							double[] region_weight = regions_weight.get(region_id).get(grp_id);

							if (val_all == null) {
								if (common_sim_header == null) {
									common_sim_header = header;
								}
								val_all = new double[region_weight.length][lines.length - 1][common_sim_header.length];

							}

							for (int r = 1; r < lines.length; r++) {
								String[] rowStr = lines[r].split(",");
								for (int c = 1; c < rowStr.length; c++) {
									int cVal = c;
									if (!common_sim_header[c].equals(header[c])) {
										System.err.printf("Warning! Column header for %s mismatch with default.\n",
												datasel_file.getAbsolutePath());
										for (int cc = 1; cc < common_sim_header.length; cc++) {
											if (common_sim_header[cc].equals(header[c])) {
												cVal = cc;
											}
										}
									}
									double col_val = Double.parseDouble(rowStr[cVal]);

									for (int ra = 0; ra < region_weight.length; ra++) {
										double[][] val = val_all[ra];
										// Time
										if (val[r - 1][0] == 0) {
											val[r - 1][0] = Double.parseDouble(rowStr[0]);
										}
										if (col_val != 0 && region_weight[ra] != 0) {
											val[r - 1][c] += region_weight[ra] * col_val;
										}
									}
								}
							}
						}
						if (flag_printProgress) {
							System.out.printf("Data extracted from %s.\n", datasel_file.getName());
						}

					} // End of for (File datasel_file : datasel_files) {

					// Print output

					for (int ra = 2; ra < regions_headers.length; ra++) {
						String region_type = regions_headers[ra].substring(2);
						double[][] val = val_all[ra - 2];

						File output_file = new File(data_dir,
								String.format(datasel_by_regions_output_fname[resPt], ent.getKey(), region_type));
						PrintWriter pWri_out = new PrintWriter(output_file);

						// Header
						pWri_out.printf("%s_%s", ent.getKey(), region_type);
						for (int i = 1; i < common_sim_header.length; i++) {
							pWri_out.print(',');
							pWri_out.print(common_sim_header[i]);
						}
						pWri_out.println();

						for (int r = 0; r < val.length; r++) {
							for (int c = 0; c < val[r].length; c++) {
								if (c == 0) {
									pWri_out.print((int) val[r][c]);
								} else {
									pWri_out.print(',');
									pWri_out.print(val[r][c]);
								}
							}
							pWri_out.println();
						}

						pWri_out.close();

						if (flag_printProgress) {
							System.out.printf("Regional specific output generated at %s\n",
									output_file.getAbsolutePath());
						}
					}

				}

				if (!keepGrpRegionCSVs) {
					File[] csvList = data_dir.listFiles(new FileFilter() {
						@Override
						public boolean accept(File pathname) {
							return pattern_datasel_csv.matcher(pathname.getName()).matches();
						}
					});

					File tar_file = new File(data_dir,
							String.format("%s.7z", String.format(datasel_by_regions_key[resPt], "ALL")));

					try {
						StaticMethods.zipFile(csvList, tar_file, true);
					} catch (IOException ex) {
						ex.printStackTrace(System.err);

					}

//					Arrays.sort(csvList, new Comparator<File>() {
//						@Override
//						public int compare(File o1, File o2) {
//							Matcher m1 = pattern_datasel_csv.matcher(o1.getName());
//							Matcher m2 = pattern_datasel_csv.matcher(o2.getName());
//							int res = 0;
//							if (m1.matches() && m2.matches()) {
//								int g = 1;
//								while (res == 0 && g < m1.groupCount() & g < m2.groupCount()) {
//									res = Integer.valueOf(m1.group(g)).compareTo(Integer.valueOf(m2.group(g)));
//									g++;
//								}
//							}
//							return res;
//						}
//					});
//
//					HashMap<Integer, ArrayList<File>> map_by_grp = new HashMap<>();
//					for (File f : csvList) {
//						Matcher m = pattern_datasel_csv.matcher(f.getName());
//						Integer grp = Integer.valueOf(m.group(1));
//						ArrayList<File> f_list = map_by_grp.get(grp);
//						if (f_list == null) {
//							f_list = new ArrayList<>();
//							map_by_grp.put(grp, f_list);
//						}
//						f_list.add(f);
//					}
//
//					for (Entry<Integer, ArrayList<File>> ent : map_by_grp.entrySet()) {
//						Integer grp = ent.getKey();
//						File[] zip_files = ent.getValue().toArray(new File[0]);
//						File tar_file = new File(data_dir, String.format("%s.7z",
//								String.format(datasel_by_regions_key[resPt]), String.format("Grp_%d", grp)));
//						StaticMethods.zipFile(zip_files, tar_file, true);
//
//					}

				}
			}

			// Group region statistics by Indigenous status and gender
			File treatment_trend_dir = new File(args[0], "Timetrend_Treatment_GrpLoc");
			File numInf_trend_dir = new File(args[0], "Timetrend_Infectious_GrpLoc");
			File adj_pop_size = new File(args[0], "Adj_pop_size.csv");
			Map<File[], Pattern> output_map;

			if (treatment_trend_dir.isDirectory()) {

				output_map = Map.ofEntries(
						Map.entry(
								new File[] {
										new File(treatment_trend_dir,
												String.format("Timetrend_Treatment_by_GrpLoc_Grp_[01234]_at_CITY.csv")),
										new File(treatment_trend_dir,
												String.format(
														"Timetrend_Treatment_by_GrpLoc_Grp_[01234]_at_REGION.csv")),
										new File(treatment_trend_dir,
												String.format(
														"Timetrend_Treatment_by_GrpLoc_Grp_[01234]_at_REMOTE.csv")), },
								Pattern.compile("Timetrend_Treatment_by_GrpLoc_Grp_([01234])_Loc_(\\d+).csv")),
						Map.entry(
								new File[] {
										new File(treatment_trend_dir,
												String.format(
														"Timetrend_Treatment_by_GrpLoc_Grp_1[01234]_at_CITY.csv")),
										new File(treatment_trend_dir,
												String.format(
														"Timetrend_Treatment_by_GrpLoc_Grp_1[01234]_at_REGION.csv")),
										new File(treatment_trend_dir,
												String.format(
														"Timetrend_Treatment_by_GrpLoc_Grp_1[01234]_at_REMOTE.csv")), },
								Pattern.compile("Timetrend_Treatment_by_GrpLoc_Grp_(1[01234])_Loc_(\\d+).csv")),
						Map.entry(
								new File[] {
										new File(treatment_trend_dir,
												String.format("Timetrend_Treatment_by_GrpLoc_Grp_[56789]_at_CITY.csv")),
										new File(treatment_trend_dir,
												String.format(
														"Timetrend_Treatment_by_GrpLoc_Grp_[56789]_at_REGION.csv")),
										new File(treatment_trend_dir,
												String.format(
														"Timetrend_Treatment_by_GrpLoc_Grp_[56789]_at_REMOTE.csv")), },
								Pattern.compile("Timetrend_Treatment_by_GrpLoc_Grp_([56789])_Loc_(\\d+).csv")),
						Map.entry(
								new File[] {
										new File(treatment_trend_dir,
												String.format(
														"Timetrend_Treatment_by_GrpLoc_Grp_1[56789]_at_CITY.csv")),
										new File(treatment_trend_dir,
												String.format(
														"Timetrend_Treatment_by_GrpLoc_Grp_1[56789]_at_REGION.csv")),
										new File(treatment_trend_dir,
												String.format(
														"Timetrend_Treatment_by_GrpLoc_Grp_1[56789]_at_REMOTE.csv")), },
								Pattern.compile("Timetrend_Treatment_by_GrpLoc_Grp_(1[56789])_Loc_(\\d+).csv")));

				System.out.println("Start analysing treatment time trend from zip...");
				Analysis_PostSim_ExtractTimeTrends.extractGrpRegionTimeTrendFromZip(
						new File(treatment_trend_dir, "Timetrend_Treatment_by_GrpLoc_ALL.csv.7z"), adj_pop_size,
						output_map);
			}

			if(numInf_trend_dir.isDirectory()) {
			output_map = Map
					.ofEntries(
							Map.entry(
									new File[] {
											new File(numInf_trend_dir,
													String.format(
															"Timetrend_NumInfect_by_GrpLoc_Grp_[01234]_at_CITY.csv")),
											new File(numInf_trend_dir,
													String.format(
															"Timetrend_NumInfect_by_GrpLoc_Grp_[01234]_at_REGION.csv")),
											new File(numInf_trend_dir, String.format(
													"Timetrend_NumInfect_by_GrpLoc_Grp_[01234]_at_REMOTE.csv")), },
									Pattern.compile("Timetrend_Infectious_by_GrpLoc_Grp_([01234])_Loc_(\\d+).csv")),
							Map.entry(
									new File[] {
											new File(numInf_trend_dir,
													String.format(
															"Timetrend_NumInfect_by_GrpLoc_Grp_1[01234]_at_CITY.csv")),
											new File(numInf_trend_dir, String.format(
													"Timetrend_NumInfect_by_GrpLoc_Grp_1[01234]_at_REGION.csv")),
											new File(numInf_trend_dir, String.format(
													"Timetrend_NumInfect_by_GrpLoc_Grp_1[01234]_at_REMOTE.csv")), },
									Pattern.compile("Timetrend_Infectious_by_GrpLoc_Grp_(1[01234])_Loc_(\\d+).csv")),
							Map.entry(
									new File[] {
											new File(numInf_trend_dir,
													String.format(
															"Timetrend_NumInfect_by_GrpLoc_Grp_[56789]_at_CITY.csv")),
											new File(
													numInf_trend_dir,
													String.format(
															"Timetrend_NumInfect_by_GrpLoc_Grp_[56789]_at_REGION.csv")),
											new File(numInf_trend_dir, String.format(
													"Timetrend_NumInfect_by_GrpLoc_Grp_[56789]_at_REMOTE.csv")), },
									Pattern.compile("Timetrend_Infectious_by_GrpLoc_Grp_([56789])_Loc_(\\d+).csv")),
							Map.entry(
									new File[] { new File(
											numInf_trend_dir,
											String.format("Timetrend_NumInfect_by_GrpLoc_Grp_1[56789]_at_CITY.csv")),
											new File(numInf_trend_dir, String.format(
													"Timetrend_NumInfect_by_GrpLoc_Grp_1[56789]_at_REGION.csv")),
											new File(numInf_trend_dir, String.format(
													"Timetrend_NumInfect_by_GrpLoc_Grp_1[56789]_at_REMOTE.csv")), },
									Pattern.compile("Timetrend_Infectious_by_GrpLoc_Grp_(1[56789])_Loc_(\\d+).csv")));

			System.out.println("Start analysing prevalence time trend from zip...");
			Analysis_PostSim_ExtractTimeTrends.extractGrpRegionTimeTrendFromZip(
					new File(numInf_trend_dir, "Timetrend_Infectious_by_GrpLoc_ALL.csv.7z"), adj_pop_size, output_map);
			}
		}
		if (flag_infection_hist_pid) {

			// Infection history for PID
			Analysis_PostSim_ExtractInfectionHistory analyseInfHist = new Analysis_PostSim_ExtractInfectionHistory(
					args);

			// PID
			int[] incl_start_grps = new int[] { 5, 6, 7, 8, 9 }; // Indigenous female
			int[] sample_time = new int[] { 6935, 7300, 7665, 8030, 8395, 8760, 9125, 9490, 9855, 10220, 10585, 10950 };

			int max_exposure = 120; // Assume won't develop PID after 4 months
			double[] event_prob_by_inf_count = new double[] { 0.14, 0.17 };
			int[] inf_count_range = new int[] { 0, 1 };

			HashMap<String, double[]> resmap = analyseInfHist.analyse(incl_start_grps, sample_time, max_exposure,
					event_prob_by_inf_count, inf_count_range);

			if (!resmap.isEmpty()) {
				File baseDir = new File(args[0]);
				File resFile = new File(baseDir, "Infection_Hist_PID.csv");
				Analysis_PostSim_ExtractInfectionHistory.generateInfectionHistCSV(resmap, sample_time, resFile);
			}

		}
	}

	public void setRegion_extract_array(String[][] region_extract_array) {
		this.region_extract_array = region_extract_array;
	}

	public void setPattern_sim(Pattern pattern_sim) {
		this.pattern_sim = pattern_sim;
	}

	public void setGrp_extract_array(String[][] grp_extract_array) {
		this.grp_extract_array = grp_extract_array;
	}

	// Index File[] corresponds to regions as specified in adj_pop_size
	public static void extractGrpRegionTimeTrendFromZip(File zipFile, File adj_pop_size,
			Map<File[], Pattern> output_pattern_map) throws IOException {

		HashMap<String, ArrayList<String[]>> map_trend_from_zip = new HashMap<>();

		Map<File, ArrayList<String>> output_builder_map = new HashMap<>();

		// Read adj pop size by region
		HashMap<Integer, HashMap<Integer, double[]>> adj_pop_size_by_region_grp = new HashMap<>();
		HashMap<Integer, HashMap<Integer, Double>> total_adj_pop_size_by_region_grp = new HashMap<>();
		String[] adj_pop_size_lines = StaticMethods.extracted_lines_from_text(adj_pop_size);

		for (int r = 1; r < adj_pop_size_lines.length; r++) {
			String[] ent = adj_pop_size_lines[r].split(",");
			Integer region_id = Integer.valueOf(ent[0]);
			Integer grp_id = Integer.valueOf(ent[1]);

			double[] adj_pop_size_by_grp = new double[ent.length - 2];
			double total_adj_pop_size = 0;

			for (int g = 2; g < ent.length; g++) {
				adj_pop_size_by_grp[g - 2] = Double.parseDouble(ent[g]);
				total_adj_pop_size += adj_pop_size_by_grp[g - 2];
			}

			HashMap<Integer, double[]> region_ent = adj_pop_size_by_region_grp.get(region_id);
			if (region_ent == null) {
				region_ent = new HashMap<>();
				adj_pop_size_by_region_grp.put(region_id, region_ent);
			}
			region_ent.put(grp_id, adj_pop_size_by_grp);

			HashMap<Integer, Double> total_region_ent = total_adj_pop_size_by_region_grp.get(region_id);
			if (total_region_ent == null) {
				total_region_ent = new HashMap<>();
				total_adj_pop_size_by_region_grp.put(region_id, total_region_ent);
			}
			total_region_ent.put(grp_id, total_adj_pop_size);
		}

		// Initialise output builder map
		for (File[] f_arr : output_pattern_map.keySet()) {
			for (File f : f_arr) {
				output_builder_map.put(f, new ArrayList<>());
			}
		}

		StaticMethods.extractedLinesFrom7Zip(zipFile, map_trend_from_zip, null);
		String[] zip_ent_names = map_trend_from_zip.keySet().toArray(new String[0]);

		for (File[] f_arr : output_pattern_map.keySet()) {
			Pattern p = output_pattern_map.get(f_arr);
			double[][][] val_all = new double[f_arr.length][][];

			for (String zip_ent_name : zip_ent_names) {
				Matcher m = p.matcher(zip_ent_name);
				if (m.matches()) {
					ArrayList<String[]> line_list = map_trend_from_zip.get(zip_ent_name);
					for (int fPt = 0; fPt < f_arr.length; fPt++) {
						File f = f_arr[fPt];
						if (f != null) {
							double[][] val = val_all[fPt];
							if (val == null) {
								val = new double[line_list.size()][line_list.get(0).length];
								// Fill in time column
								for (int r = 1; r < line_list.size(); r++) {
									val[r][0] = Integer.parseInt(line_list.get(r)[0]);
								}
								// Initialise output builder header
								if (output_builder_map.get(f).isEmpty()) {
									StringBuilder header = new StringBuilder();
									header.append("Time");
									for (int c = 1; c < line_list.get(0).length; c++) {
										header.append(',');
										header.append(line_list.get(0)[c]);
									}
									output_builder_map.get(f).add(header.toString());
								}
								val_all[fPt] = val;
							}
						}
					}

					Integer grp_id = Integer.valueOf(m.group(1));
					Integer region_id = Integer.valueOf(m.group(2));
					double[] adj_pop_size_by_grp = adj_pop_size_by_region_grp.get(region_id).get(grp_id);
					double total_adj_pop_size = total_adj_pop_size_by_region_grp.get(region_id).get(grp_id);

					for (int r = 1; r < line_list.size(); r++) {
						String[] line = line_list.get(r);
						for (int c = 1; c < line.length; c++) {
							int baseVal = Integer.parseInt(line[c]);
							for (int fPt = 0; fPt < f_arr.length; fPt++) {
								if (f_arr[fPt] != null) {
									double[][] val = val_all[fPt];
									if (val != null) {
										val[r][c] += baseVal * adj_pop_size_by_grp[fPt] / total_adj_pop_size;
									}
								}
							}
						}
					}
				}
			}

			for (int fPt = 0; fPt < f_arr.length; fPt++) {
				File f = f_arr[fPt];
				double[][] val = val_all[fPt];

				// Fill in output builder
				if (val != null) {
					for (int r = 1; r < val.length; r++) {
						StringBuilder lineBuilder = new StringBuilder();
						lineBuilder.append((int) val[r][0]);
						for (int c = 1; c < val[r].length; c++) {
							lineBuilder.append(',');
							lineBuilder.append(val[r][c]);
						}
						output_builder_map.get(f).add(lineBuilder.toString());
					}

					// Printing of outputs
					PrintWriter pWri_out = new PrintWriter(f);
					for (String line : output_builder_map.get(f)) {
						pWri_out.println(line);
					}
					pWri_out.close();

					System.out.printf("Time trend extracted to %s\n", f.getAbsolutePath());

				}
			}

		}

	}
}
