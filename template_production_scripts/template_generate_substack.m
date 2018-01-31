% generate renderer collection struct 

rcfrom.stack          = 'v12_acquire_LC_merged';
rcfrom.owner          ='flyTEM';
rcfrom.project        = 'FAFB00';
rcfrom.service_host   = '10.40.3.162:8080';
rcfrom.baseURL        = ['http://' rcfrom.service_host '/render-ws/v1'];
rcfrom.verbose        = 0;


rcto.stack          = 'FAFB_section_one_acquire';
rcto.owner          ='flyTEM';
rcto.project        = 'test_solver';
rcto.service_host   = '10.40.3.162:8080';
rcto.baseURL        = ['http://' rcto.service_host '/render-ws/v1'];
rcto.verbose        = 0;



[err] = copy_renderer_section(1, rcfrom, rcto, '/scratch/khairyk/');

set_renderer_stack_state_complete(rcto);

%%
rcto.stack          = 'FAFB_three_section_slab_acquire';
rcto.owner          ='flyTEM';
rcto.project        = 'test_solver';
rcto.service_host   = '10.40.3.162:8080';
rcto.baseURL        = ['http://' rcto.service_host '/render-ws/v1'];
rcto.verbose        = 0;



[err] = copy_renderer_section([1:3], rcfrom, rcto, '/scratch/khairyk/');

set_renderer_stack_state_complete(rcto);
%%
% generate renderer collection struct 

rcfrom.stack          = 'v13_montage';
rcfrom.owner          ='flyTEM';
rcfrom.project        = 'FAFB00';
rcfrom.service_host   = '10.40.3.162:8080';
rcfrom.baseURL        = ['http://' rcfrom.service_host '/render-ws/v1'];
rcfrom.verbose        = 0;


rcto.stack          = 'FAFB_montage_three_section';
rcto.owner          ='flyTEM';
rcto.project        = 'test_solver';
rcto.service_host   = '10.40.3.162:8080';
rcto.baseURL        = ['http://' rcto.service_host '/render-ws/v1'];
rcto.verbose        = 0;



[err] = copy_renderer_section([1:3], rcfrom, rcto, '/scratch/khairyk/');

set_renderer_stack_state_complete(rcto);