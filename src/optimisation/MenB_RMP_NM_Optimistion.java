package optimisation;

import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.util.Arrays;
import java.util.Comparator;
import java.util.HashMap;
import java.util.Properties;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;

import org.apache.commons.math3.optim.InitialGuess;
import org.apache.commons.math3.optim.MaxEval;
import org.apache.commons.math3.optim.PointValuePair;
import org.apache.commons.math3.optim.nonlinear.scalar.GoalType;
import org.apache.commons.math3.optim.nonlinear.scalar.MultivariateFunctionMappingAdapter;
import org.apache.commons.math3.optim.nonlinear.scalar.ObjectiveFunction;
import org.apache.commons.math3.optim.nonlinear.scalar.noderiv.NelderMeadSimplex;
import org.apache.commons.math3.optim.nonlinear.scalar.noderiv.SimplexOptimizer;

public class MenB_RMP_NM_Optimistion {

	private final HashMap<String, String> cross_ref_sample_range;
	private final String[] param_to_opt;

	private final int[] opt_time_range;
	private final String[] opt_outcome_csv;
	private final double[][] opt_setting;

	private final String path_dirName;
	private final String path_region_path;
	private final String path_grp_size;

	private final String path_seed_dir;
	private final String[] seed_file_lines;
	private final String[] seed_file_header;
	
	private final double[][] param_boundaries;
	
	
	private final double RELATIVE_TOLERANCE = 1e-5;
	private final double ABSOLUTE_TOLERANCE = 1e-10;
	private final MaxEval maxVal = MaxEval.unlimited(); // new MaxEval(numEval);

	private boolean printProgress = false;
	
	

	public MenB_RMP_NM_Optimistion(String dirName) throws IOException {
		this.path_dirName = dirName;

		File file_opt_setting = new File(dirName, "optSetting.prop");
		FileInputStream fIS = new FileInputStream(file_opt_setting);
		Properties prop = new Properties();
		prop.loadFromXML(fIS);
		fIS.close();

		String param_to_opt_str = prop.getProperty("PROP_PARAM_TO_OPT");

		param_to_opt = param_to_opt_str.split(",");
		HashMap<String, double[]> default_sample_range = new HashMap<>();
		cross_ref_sample_range = new HashMap<>();

		for (String param : param_to_opt) {
			String ent = prop.getProperty(String.format("PROP_PARAM_SETTING_%s", param));
			if (ent != null) {
				ent.replaceAll("\\s", "");
				String[] sp = ent.split(",");
				default_sample_range.put(param, new double[] { Double.parseDouble(sp[0]), Double.parseDouble(sp[1]) });
				if (sp.length > 2) {
					cross_ref_sample_range.put(param, sp[2]);
				}
			}
		}

		path_seed_dir = prop.getProperty("PPOP_PATH_SEED_DIR");
		path_region_path = prop.getProperty("PROP_PATH_REGION_MAPPING");
		path_grp_size = prop.getProperty("PROP_PATH_GRP_SIZE");

		opt_time_range = (int[]) util.PropValUtils.propStrToObject(prop.getProperty("PROP_OPT_TIME_RANGE"),
				int[].class);
		opt_outcome_csv = prop.getProperty("PROP_OPT_OUTCOME_CSV").replaceAll("\\s", "").split(",");

		opt_setting = (double[][]) util.PropValUtils.propStrToObject(prop.getProperty("PROP_OPT_SETTING"),
				double[][].class);

		// Set up initial parameter array
		seed_file_lines = util.Util_7Z_CSV_Entry_Extract_Callable.extracted_lines_from_text(
				new File(new File(dirName, path_seed_dir), String.format("%s.csv", path_seed_dir)));
		seed_file_header = seed_file_lines[0].split(",");
		
		
		param_boundaries = new double[2][param_to_opt.length];
		for (int i = 0; i < param_to_opt.length; i++) {
			double[] range = default_sample_range.get(param_to_opt[i]);
			param_boundaries[0][i] = range[0];
			param_boundaries[1][i] = range[1];
		}
		

	}

	public void runOptimisation() {
		ExecutorService exec = Executors.newFixedThreadPool(seed_file_lines.length-1);			
		for (int p = 1; p < seed_file_lines.length; p++) {
			String[] seed_file_def_val = seed_file_lines[p].split(",");
			HashMap<String, Double> init_value = new HashMap<>();
			for (int i = 0; i < seed_file_header.length; i++) {
				init_value.put(seed_file_header[i], Double.valueOf(seed_file_def_val[i]));
			}			
			double[] param_init = new double[param_to_opt.length];
			for (int i = 0; i < param_to_opt.length; i++) {				
				param_init[i] = init_value.get(param_to_opt[i]).doubleValue();

			}
			// Adjust for cross reference
			for (int i = 0; i < param_to_opt.length; i++) {
				if (cross_ref_sample_range.containsKey(param_to_opt[i])) {
					param_init[i] = param_init[i] / init_value.get(cross_ref_sample_range.get(param_to_opt[i]));
				}
			}			
			String[] path = new String[] { path_dirName,String.format("%s_%d", path_seed_dir, p-1),	path_region_path,path_grp_size};
			
			// Objective function
			ResidualFunc_RMP func = new ResidualFunc_RMP(path, new String[][] { seed_file_header, seed_file_def_val },
					param_to_opt, cross_ref_sample_range, opt_outcome_csv, opt_setting, opt_time_range);
			func.setPrintProgess(printProgress);
			
			// Set up simplex			
			
			MultivariateFunctionMappingAdapter wrapper = new MultivariateFunctionMappingAdapter(func, param_boundaries[0],
					param_boundaries[1]);
			
			
			ObjectiveFunction objFunc = new ObjectiveFunction(wrapper);
			
			Runnable runnable_opt = new Runnable() {				
				@Override
				public void run() {
					InitialGuess initial_guess;
					final NelderMeadSimplex simplex;

					initial_guess = new InitialGuess(wrapper.boundedToUnbounded(param_init));
					simplex = new NelderMeadSimplex(param_init.length);

					SimplexOptimizer optimizer = new SimplexOptimizer(RELATIVE_TOLERANCE, ABSOLUTE_TOLERANCE);
					
					try {
						PointValuePair pV;
						System.out.printf("Optimisation start Max Eval = %d.\n", maxVal.getMaxEval());

						pV = optimizer.optimize(objFunc, simplex, GoalType.MINIMIZE, initial_guess, maxVal);
						double[] point = wrapper.unboundedToBounded(pV.getPoint());

						StringBuilder pt_str = new StringBuilder();
						for (double pt : point) {
							if (pt_str.length() != 0) {
								pt_str.append(',');
							}
							pt_str.append(String.format("%.5f", pt));
						}

						System.out.printf("Optimisation Completed.\nP = [%s], V = %f\n", pt_str.toString(), pV.getValue());

					} catch (org.apache.commons.math3.exception.TooManyEvaluationsException ex) {
						System.out.printf("Eval limit of (ex.getMax=%d) reached.\nSimplex (bounded):\n", ex.getMax());

						PointValuePair[] res = simplex.getPoints();
						Arrays.sort(res, new Comparator<PointValuePair>() {
							@Override
							public int compare(PointValuePair o1, PointValuePair o2) {
								return Double.compare(o1.getValue(), o2.getValue());
							}

						});

						for (PointValuePair pV : res) {
							double[] point = wrapper.unboundedToBounded(pV.getPoint());

							StringBuilder pt_str = new StringBuilder();
							for (double pt : point) {
								if (pt_str.length() != 0) {
									pt_str.append(',');
								}
								pt_str.append(String.format("%.5f", pt));
							}

							System.out.printf("P = [%s], V = %f\n", pt_str.toString(), pV.getValue());

						}

					}
					
					
				}
			};
			
			if(seed_file_lines.length == 2) {
				runnable_opt.run();
			}else {
				exec.execute(runnable_opt);
			}

		}
		
		if(seed_file_lines.length > 2) {
			exec.shutdown();
			try {
				if (!exec.awaitTermination(2, TimeUnit.DAYS)) {
					System.err.println("Thread time-out!");
				}
			} catch (InterruptedException e) {				
				e.printStackTrace(System.err);
			}
		}

	}
	
	public void setPrintProgress(boolean printProgress) {
		this.printProgress = printProgress;
	}

}
