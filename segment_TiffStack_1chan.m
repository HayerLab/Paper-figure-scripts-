root = 'Z:\Arnold\230705_TIE'; 
output_path = ([root, filesep, 'masks']);


files = getFilenames(root, '.tif');

for i = 1:size(files, 1)

    stack=double(readTIFFstack([root, filesep, files{i,1}]));
     n_frames = size(stack, 3);

    maskFinal = {}; 
        cellCoors = {}; 
        outlined_stack = {}; 
        masked_stack = cell(1, n_frames);

   

        for frame = 1:n_frames
            image = stack(:,:,frame);
            [mask, coordinates] = getCellMaskCyto_3(image, 0.25, 3,0,0);
            
            maskFinal{frame} = mask;
            cellCoors{frame} = coordinates;
           
            
            masked_image = double(image);
            masked_image(~mask) = nan; 
            
            masked_stack{frame} = uint16(masked_image);
            
           
            outline = bwperim(mask, 8);
        %     colored_outline = ?
            img_outline{frame} = outline;
            
            outlined_stack{frame} = DrawMaskOutline(stack(:,:,frame), mask);
            
             imwrite(img_outline{frame}, [output_path, files{i,1}, '_outline.tif'], 'WriteMode','append','Compression','none');
            imwrite(masked_stack{frame}, [output_path, files{i,1}, '_masked_img.tif'], 'WriteMode','append','Compression','none');
%             imwrite(outlined_stack{frame}, [output_path, num2str(char(file_name(file))), '-0', num2str(site), '_outlined_img_sf45_pj6.tif'], 'WriteMode','append','Compression','none');
       
        
    
               
        save([output_path, files{i,1}, '_mask.mat'],'maskFinal','cellCoors','outlined_stack','-v7.3'); 
%         save([output_path, num2str(char(file_name(file))),  '-0', num2str(site), '_bleach.mat'],'image_mean', 'masked_image_mean','-v7.3'); 
    end
end

   



