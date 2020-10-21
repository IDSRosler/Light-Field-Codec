DATASET_DIR = 'C:\Users\Cristian\Documents\UFPel\lfcodec\Datasets';
RESULTS_DIR = 'C:\Users\Cristian\Documents\UFPel\lfcodec\Results';
DATASETS = [
     "Bikes",
    "Danger_de_Mort",   
    "Fountain_Vincent2",
    "Stone_Pillars_Outside",
    "poznanlab1_TAU",
    "tarot",
];

TRANSFORMS = ["DCT" "DST"];
QUANTIZATIONS = ["0.1" "0.25" "0.5" "1" "2" "5" "7" "10" "15" "20" "50" "100"];
REF_FOLDER_FMT = "%s\\%s\\";
REC_FOLDER_FMT = '%s\\%s_15_15_13_13_%s_%s_%s_%s_%s\\';
OUT_FILE_FMT = '%s_15_15_13_13_%s_%s_%s_%s_%s.csv';

TOTAL = length(DATASETS) * length(TRANSFORMS) * length(QUANTIZATIONS);
run = 0;
for i_dataset=1:length(DATASETS)
    dataset = DATASETS(i_dataset);
    for i_transform=1:length(TRANSFORMS)
        tx = TRANSFORMS(i_transform);
        for i_quantization=1:length(QUANTIZATIONS)
            qt = QUANTIZATIONS(i_quantization);
            ref_folder = sprintf(REF_FOLDER_FMT, DATASET_DIR, dataset);
            rec_folder = sprintf(REC_FOLDER_FMT, ...
                                 RESULTS_DIR, dataset, qt,qt,qt,qt,tx);
            out_file = sprintf(OUT_FILE_FMT, dataset, qt,qt,qt,qt,tx);
            percentage = run / TOTAL;
            try
                fprintf("[%5.2f%%] Computing %s ... ", percentage, out_file);
                compute(ref_folder, rec_folder, out_file);
                fprintf("OK!\n");
            catch
                fprintf("Something wrong happened :(\n");
            end
            run = run + 1;
        end
    end
end

