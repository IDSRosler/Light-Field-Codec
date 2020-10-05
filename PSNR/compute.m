function compute(folder_orig, folder_modf, outputFile)

    files = dir(strcat(folder_modf,'*.ppm'));

    fileID = fopen(outputFile, 'w');
    fprintf(fileID,'View, Y_PSNR, YUV_PSNR, Y_SSIM\n');

    for it=1:length(files)
        view_name = files(it).name;

        file_orig = strcat(folder_orig, view_name);
        file_modf = strcat(folder_modf, view_name);

        im_orig = imread(file_orig);
        im_modif = imread(file_modf);

        ref = bitshift(im_orig,-6); rec = bitshift(im_modif,-6);

        [Y_PSNR, YUV_PSNR, Y_SSIM] = QM(ref, rec, 10, 10);
        
        % result = [view_name(1:end-4), Y_PSNR, YUV_PSNR, Y_SSIM];
        % disp([view_name(1:end-4), ' ', num2str(Y_PSNR), ' ', num2str(YUV_PSNR), ' ', num2str(Y_SSIM)]);

        fprintf(fileID, '%s, %f, %f, %f\n', view_name(1:end-4), Y_PSNR, YUV_PSNR, Y_SSIM);
    end
    
    fclose(fileID);
end