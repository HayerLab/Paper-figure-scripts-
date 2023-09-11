function [] = segment_TIFFstack(filepath, mincellsize)
stack=double(readTIFFstack([filepath]));
s=1; 
frameNum=1;
for i =1:size(stack,3)
    
    image = stack(:,:,1); 
    
    
    
   
    
    
    
    
  frameNum=frameNum+1;   
end 



end

for file = files_to_analyze
    masked_mean_bleaching = [];
    bleach_mean = [];
    
    for site = sites(file, 1:n_sites_per_file(file))
        disp([num2str(char(file_name(file))), '-0', num2str(site)])
        stack_path = [TIFFStacks_path, num2str(char(file_name(file))), '-0', num2str(site), '.tif'];
        warning off
        stack = readTIFFstack(stack_path);

        n_frames = size(stack, 3);
        maskFinal = {};
        cellCoors = {};
        im_actin = {};
        
        outlined_stack = {};
        masked_stack = cell(1, n_frames);
        masked_image_mean = [];
        image_mean = [];

        for frame = 1:n_frames
            image = stack(:,:,frame);
            [mask, coordinates] = getCellMaskCyto_3(image);
            
            maskFinal{frame} = mask;
            cellCoors{frame} = coordinates;
            im_actin{frame} = image;
            
            masked_image = double(image);
            masked_image(~mask) = nan; 
            
            masked_stack{frame} = uint16(masked_image);
            
            % For bleaching estimation:
            masked_image_mean(frame) = nanmean(vect(masked_stack{frame}));
            image_mean(frame) = nanmean(vect(stack(:,:,frame)));
            
            outline = bwperim(mask, 8);
        %     colored_outline = ?
            img_outline{frame} = outline;
            
            outlined_stack{frame} = DrawMaskOutline(stack(:,:,frame), mask);
            
%             imwrite(img_outline{frame}, [output_path, num2str(char(file_name(file))), '-0', num2str(site), '_outline.tif'], 'WriteMode','append','Compression','none');
%             imwrite(masked_image{frame}, [output_path, num2str(char(file_name(file))), '-0', num2str(site), '_masked_img.tif'], 'WriteMode','append','Compression','none');
%             imwrite(outlined_stack{frame}, [output_path, num2str(char(file_name(file))), '-0', num2str(site), '_outlined_img_sf45_pj6.tif'], 'WriteMode','append','Compression','none');
        end
        
        masked_mean_bleaching(:, site, file) = masked_image_mean / masked_image_mean(1);
%         plot(1:n_frames, masked_mean_bleaching(1:n_frames, site, file), ':'); hold on    
               
        save([masks_path, num2str(char(file_name(file))),  '-0', num2str(site), '_mask.mat'],'maskFinal','cellCoors','im_actin','outlined_stack','-v7.3'); 
%         save([output_path, num2str(char(file_name(file))),  '-0', num2str(site), '_bleach.mat'],'image_mean', 'masked_image_mean','-v7.3'); 
    end
end