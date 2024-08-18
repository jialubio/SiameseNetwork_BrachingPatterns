function [plaintext_binary, plaintext] = RandPlainText_Generator(plaintext_length,char_list,binary_list)
% plaintext_length = 3; 
rand_vec = rand(plaintext_length,1);
plaintext_char = ceil(36*rand_vec);
plaintext_cell = char_list(plaintext_char);
plaintext_binary = binary_list(plaintext_char,:);

for i = 1:plaintext_length,
    if i ==1 ,
        plaintext = plaintext_cell{i};
    else
        plaintext = strcat(plaintext,plaintext_cell{i});
    end
end
end

