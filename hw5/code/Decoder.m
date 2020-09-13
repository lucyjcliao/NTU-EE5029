function u_decode = Decoder(q,trun_len)
% [Viterbi decoder]
% Argument:
%   - q: receive signals
%   - trun_len: truncation length
% Return:
%   - u_decode: decoded codeword
u_decode = [];
state_buffer = [0 0 0 0];                                                   % Initial state
cost = zeros(16,trun_len);
path = Inf(16,trun_len);
u_decode = [];
i = 0;
q = [q zeros(1,trun_len*2)];

while length(u_decode) < 1000
    i = i + 1;
    r_code = [q(2*(i-1)+1) q(2*(i-1)+2)];                                   % receive codeword
    next_state_buffer = [];
    new_path = Inf(16,trun_len);                                            % create a new path buffer
    cost = [Inf(16,1) cost(:,1:end-1)];                                     % create a new column buffer
    
    for j = 1:size(state_buffer,1)
        current_state = state_buffer(j,:);

        % input = 0
        [codeword_0,next_state_0] = ConvolutionalCode(0,current_state);
        mol_0 = Modulation(codeword_0);
        cost_0 = sum((r_code-mol_0).^2);
        if cost(bi2de(next_state_0,'left-msb')+1,1)==Inf || cost(bi2de(next_state_0,'left-msb')+1,1) >= cost(bi2de(current_state,'left-msb')+1,2)+cost_0
            cost(bi2de(next_state_0,'left-msb')+1,1) = cost(bi2de(current_state,'left-msb')+1,2)+cost_0;
            new_path(bi2de(next_state_0,'left-msb')+1,:) = [0,path(bi2de(current_state,'left-msb')+1,1:end-1)];
        end
        next_state_buffer = [next_state_buffer;next_state_0];

        % input = 1
        [codeword_1,next_state_1] = ConvolutionalCode(1,current_state);
        mol_1 = Modulation(codeword_1);
        cost_1 = sum((r_code-mol_1).^2);
        if cost(bi2de(next_state_1,'left-msb')+1,1)==Inf || cost(bi2de(next_state_1,'left-msb')+1,1) >= cost(bi2de(current_state,'left-msb')+1,2)+cost_1
            cost(bi2de(next_state_1,'left-msb')+1,1) = cost(bi2de(current_state,'left-msb')+1,2)+cost_1;
            new_path(bi2de(next_state_1,'left-msb')+1,:) = [1,path(bi2de(current_state,'left-msb')+1,1:end-1)];
        end
        next_state_buffer = [next_state_buffer;next_state_1];

    end
    
    % when depth > 12
    if sum(cost(:,end)) ~= 0
        best = find(cost(:,1)==min(cost(:,1)));
        u_decode(end+1) = mode(new_path(best,end));
    end
    
    % update data
    path = new_path;
    state_buffer = unique(next_state_buffer,'rows');
    
end
end 
