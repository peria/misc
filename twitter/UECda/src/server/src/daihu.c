#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include<time.h>
#include<string.h>
#include<strings.h>
#include<unistd.h>
#include<sys/socket.h>
#include<arpa/inet.h>
#include<sys/param.h>
#include<string.h>
#include<sys/fcntl.h>
#include<sys/types.h>
#include<sys/socket.h>
#include<netinet/in.h>
#include<netdb.h>
#include<pthread.h>
#include<errno.h>
#include<math.h>
#include<time.h>

#include<X11/Intrinsic.h>
#include<X11/StringDefs.h>
#include<X11/Xaw/Label.h>
#include<X11/xpm.h>

#include "xdaihu.h"
#include "external.h"
#include "tn_protocol.h"

#include "statistics.h"

// #include "config.h"

int main(int argc, char *argv[]){

	/* setting of valuances */

	int i,j,k;		      // roop
	char *dat1;		      // roop
	char *dat2;		      // roop
	char dat3[256];		      // roop
	char dat4[256];		      // roop

		/* controle game */
	int flag_wait_type=3;
	int graph_mode=0;
	int g_flag=0;
	int tmp;

		/* for initialize_windows2 */
	int accept_flag=0;

		/* for cards */
	int work_card[8][15]={0};       // i.e. submitted card etc...
	int stage_card[8][15]={0};       // i.e. submitted card etc...
	int old_stage_card[8][15]={0};       // i.e. submitted card etc...
	int players_card[5][8][15]={0}; // players_card
	int tmp_card[2][8][15]={0}; // players_card
	int number_of_card;	     // number_of_card temporal
	int status_of_submitted_card[4]={0}; // now card
	int status_of_stages_card[4]={-1};  // cards on the stage

	int sekijun[5]; // sekijun
	int mibun[5];  // cards on the stage
	int human[5]={0,0,0,0,0};  // cards on the stage
	unsigned int point[5][15000]={0,0,0,0,0};  // cards on the stage
	int error;		      // error check
	
	int now_number_of_games;	// save number of played games

		/* flags */	
	int now_pass[5]={0,0,0,0,0};
	int flash=0;   // flash
	int last_player; 
	int now_muki;
	int now_player;  // now playing player number
	int now_kakumei; // flag of kakumei
	int now_11back;  // flag of 11 back
	int now_number_of_goal=0;  // number of person whose card is none
	int goal_flag[5]={0};      // flags of person  
	int now_jocker_single=0;   // flag of jocker
	int now_shibari=0;   // flag of jocker
	int number_renzoku_pass=0;       // flags of passed person

		/* for tcp */
	int sockfd, client_sockfd[5];
	int port_number=42485;
	int protocol_version=20070;
	struct sockaddr_in wait_addr; // waiting port
	struct sockaddr_in client_addr[5]; // port for each clients
	int client_len[5]; // waiting port
	fd_set target_fds;
	fd_set org_target_fds;
	struct timeval waitval;

		/* for log */
	FILE *fp,*fp2;
	char cfg_file[100]="tndhms.cfg";
	int debug=0;     // debug_flag

		/* for window */
	Display *display;
	int root;
	int screen;
	GC gc;
	unsigned long foreground;
	unsigned long background;
	XSetWindowAttributes attr;
	int err;

	char XPM_PATH[80]="\0";
	char XPM_CARD[100];
	char XPM_TEFUDA[100];
	char XPM_TEFUDA2[100];
	char XPM_MIBUN[100];

	c_window_info win_info;
	g_window_info g_win_info;
	int WINDOW_TYPE=0;
	int GRAPH_WINDOW=0;
	char player_name[5][15];

		/* for rule */
	int RULE_KAKUMEI=1;
	int RULE_SHIBARI=1;
	int RULE_KINSOKU=0;
	int RULE_CHANGE=1;
	int RULE_KAIDAN=1;
	int RULE_5TOBI=0;
	int RULE_6REVERS=0;
	int RULE_8GIRI=1;
	int RULE_11BACK=0;
	int RULE_SEKIGAE=1;
	int RULE_SEKIGAE_NUM=3;
	int RAND_TYPE=0;
	int GAME_NUMBER=100;
	int FLASH_MIBUN_NUMBER=GAME_NUMBER;
	
	/* for player statistics */
	int count_turn;
	struct playerStatistics ps[5];
	for(i=0;i<=4;i++){
		bzero(&ps[i], sizeof(struct playerStatistics));
	}
	int sum_of_turn = 0;
	int game_count = 0;

	/************************************/
	/*  setting initial value	    */
	/************************************/

	for(i=0;i<=4;i++){
		sprintf(player_name[i],"Player %i\0",i+1);
	}

	/************************************/
	/*  setting for argument part1	    */
	/************************************/

	for(i=1;i<=argc-1;i++){
		if(strcmp(argv[i],"-v")==0){ /* print version */
			printf("tndhms version 0.29\n");
			return 0;
		}else if((strcmp(argv[i],"-h")==0)){ /* print help message */
			printf("tndhms [-vh] [-p port_number] \n");
			printf("   -v version\n");
			printf("   -h help\n");
			printf("   -c config_file\n");
			printf("   -p port_number\n");
			return 0;
		}else if((strcmp(argv[i],"-c")==0)){ /* config file name */
			if(i+1<=argc-1){
				fp=fopen(argv[i+1],"r");
				if((strlen(argv[i+1])>=90) || fp==NULL){
					printf("Bad file name. \n");
					return 0;	
				}else{
					strcpy(cfg_file,argv[i+1]);
					printf("a use config_file is %s\n", cfg_file);
				}
				fclose(fp);
			}else{
				printf("Bad file name. \n");
				return 0;	
			}
		}
	}

	/************************************/
	/*  read config file	            */
	/************************************/
	fp=fopen(cfg_file,"r");
	if(fp==NULL){
		printf("I can open no config file.\n");
		return 0;
	};	
	error=0;
	k=0;
	while(fgets(dat3,256,fp)!= NULL){
		k++;
		if(dat3[0]=='#'){
			strcpy(dat3," \0");
		}
		strcpy(dat4,dat3);
		dat2=strtok(dat3,"#");
		strupr2(dat2);
		dat1=strtok(dat2," \r\t\n");
		if(dat1!=NULL){
			if(strcmp(dat1,"WINDOW_TYPE")==0){
				dat1=strtok(NULL," \r\t\n");
				if(strcmp(dat1,"SMALL")==0){
					WINDOW_TYPE=0;	
				}else if(strcmp(dat1,"BIG")==0){
					WINDOW_TYPE=1;	
				}else if(strcmp(dat1,"CONSOLE")==0){
					WINDOW_TYPE=2;	
				}else if(strcmp(dat1,"MIDDLE")==0){
					WINDOW_TYPE=4;	
				}else if(strcmp(dat1,"MIDDLE_EXTEND")==0){
					WINDOW_TYPE=3;	
				}else if(strcmp(dat1,"MIDDLE_XGA")==0){
					WINDOW_TYPE=5;	
				}else if(strcmp(dat1,"NEW_NORMAL")==0){
					WINDOW_TYPE=6;	
				}else{
					error=1;
				}
			}else if(strcmp(dat1,"GRAPH_WINDOW")==0){
				dat1=strtok(NULL," \r\t\n");
				if(strcmp(dat1,"NONE")==0 || strcmp(dat1,"NO")==0){
					GRAPH_WINDOW=0;	
				}else if(strcmp(dat1,"BIG")==0){
					GRAPH_WINDOW=1;	
				}else if(strcmp(dat1,"MIDDLE")==0){
					GRAPH_WINDOW=2;	
				}else{
					error=1;
				}
			}else if(strcmp(dat1,"GAME_NUMBER")==0){
				dat1=strtok(NULL," \r\t\n");
				if(isint(dat1)){
					GAME_NUMBER=atoi(dat1);	
				}else{
					error=1;
				}
			}else if(strcmp(dat1,"FLASH_MIBUN_NUMBER")==0){
				dat1=strtok(NULL," \r\t\n");
				if(isint(dat1)){
					FLASH_MIBUN_NUMBER=atoi(dat1);	
				}else{
					error=1;
				}
			}else if(strcmp(dat1,"PORT_NUMBER")==0){
				dat1=strtok(NULL," \r\t\n");
				if(isint(dat1)){
					port_number=atoi(dat1);	
				}else{
					error=1;
				}
			}else if(strcmp(dat1,"RULE_KAKUMEI")==0){
				dat1=strtok(NULL," \r\t\n");
				if(strcmp(dat1,"NO")==0){
					RULE_KAKUMEI=0;	
				}else if(strcmp(dat1,"YES")==0){
					RULE_KAKUMEI=1;	
				}else{
					error=1;
				}
			}else if(strcmp(dat1,"RULE_SHIBARI")==0){
				dat1=strtok(NULL," \r\t\n");
				if(strcmp(dat1,"NO")==0){
					RULE_SHIBARI=0;	
				}else if(strcmp(dat1,"YES")==0){
					RULE_SHIBARI=1;	
				}else{
					error=1;
				}
			}else if(strcmp(dat1,"RULE_KINSOKU")==0){
				dat1=strtok(NULL," \r\t\n");
				if(strcmp(dat1,"NO")==0){
					RULE_KINSOKU=0;	
				}else if(strcmp(dat1,"YES")==0){
					RULE_KINSOKU=1;	
				}else{
					error=1;
				}
			}else if(strcmp(dat1,"RULE_KAIDAN")==0){
				dat1=strtok(NULL," \r\t\n");
				if(strcmp(dat1,"NO")==0){
					RULE_KAIDAN=0;	
				}else if(strcmp(dat1,"YES")==0){
					RULE_KAIDAN=1;	
				}else{
					error=1;
				}
			}else if(strcmp(dat1,"RULE_5TOBI")==0){
				dat1=strtok(NULL," \r\t\n");
				if(strcmp(dat1,"NO")==0){
					RULE_5TOBI=0;	
				}else if(strcmp(dat1,"YES")==0){
					RULE_5TOBI=1;	
				}else{
					error=1;
				}
			}else if(strcmp(dat1,"RULE_6REVERS")==0){
				dat1=strtok(NULL," \r\t\n");
				if(strcmp(dat1,"NO")==0){
					RULE_6REVERS=0;	
				}else if(strcmp(dat1,"YES")==0){
					RULE_6REVERS=1;	
				}else{
					error=1;
				}
			}else if(strcmp(dat1,"RULE_8GIRI")==0){
				dat1=strtok(NULL," \r\t\n");
				if(strcmp(dat1,"NO")==0){
					RULE_8GIRI=0;	
				}else if(strcmp(dat1,"YES")==0){
					RULE_8GIRI=1;	
				}else{
					error=1;
				}
			}else if(strcmp(dat1,"RULE_11BACK")==0){
				dat1=strtok(NULL," \r\t\n");
				if(strcmp(dat1,"NO")==0){
					RULE_11BACK=0;	
				}else if(strcmp(dat1,"YES")==0){
					RULE_11BACK=1;	
				}else{
					error=1;
				}
			}else if(strcmp(dat1,"RULE_CHANGE")==0){
				dat1=strtok(NULL," \r\t\n");
				if(strcmp(dat1,"NO")==0){
					RULE_CHANGE=0;	
				}else if(strcmp(dat1,"YES")==0){
					RULE_CHANGE=1;	
				}else{
					error=1;
				}
			}else if(strcmp(dat1,"RULE_SEKIGAE")==0){
				dat1=strtok(NULL," \r\t\n");
				if(strcmp(dat1,"NO")==0){
					RULE_SEKIGAE=0;	
				}else if(strcmp(dat1,"YES")==0){
					RULE_SEKIGAE=1;	
				}else{
					error=1;
				}
			}else if(strcmp(dat1,"RULE_SEKIGAE_NUM")==0){
				dat1=strtok(NULL," \r\t\n");
				if(isint(dat1)){
					RULE_SEKIGAE_NUM=atoi(dat1);	
				}else{
					error=1;
				}
			}else if(strcmp(dat1,"PLAYER1_NAME")==0){
				dat1=strtok(dat4," \r\t\n");
				dat1=strtok(NULL," \r\t\n");
				if(dat1!=NULL){
					snprintf(player_name[0],9,"%s",dat1);
				}
			}else if(strcmp(dat1,"PLAYER2_NAME")==0){
				dat1=strtok(dat4," \r\t\n");
				dat1=strtok(NULL," \r\t\n");
				if(dat1!=NULL){
					snprintf(player_name[1],9,"%s",dat1);
				}
			}else if(strcmp(dat1,"PLAYER3_NAME")==0){
				dat1=strtok(dat4," \r\t\n");
				dat1=strtok(NULL," \r\t\n");
				if(dat1!=NULL){
					snprintf(player_name[2],9,"%s",dat1);
				}
			}else if(strcmp(dat1,"PLAYER4_NAME")==0){
				dat1=strtok(dat4," \r\t\n");
				dat1=strtok(NULL," \r\t\n");
				if(dat1!=NULL){
					snprintf(player_name[3],9,"%s",dat1);
				}
			}else if(strcmp(dat1,"PLAYER5_NAME")==0){
				dat1=strtok(dat4," \r\t\n");
				dat1=strtok(NULL," \r\t\n");
				if(dat1!=NULL){
					snprintf(player_name[4],9,"%s",dat1);
				}
			}else if(strcmp(dat1,"XPM_DIR_PATH")==0){
				dat1=strtok(dat4," \r\t\n");
				dat1=strtok(NULL," \r\t\n");
				if(dat1!=NULL){
					snprintf(XPM_PATH,100,"%s/",dat1);
				}
			}else if(strcmp(dat1,"RAND_TYPE")==0){
				dat1=strtok(dat4," \r\t\n");
				dat1=strtok(NULL," \r\t\n");
				if(strcmp(dat1,"ISO_NORMAL")==0){
					RAND_TYPE=0;	
				}else if(strcmp(dat1,"MT")==0){
					RAND_TYPE=1;	
				}
			}

			if(error){
				printf("I found error in line %i\n",k);
				error=0;
			}
		}
	}
	fclose(fp);

	/************************************/
	/*  setting for argument part2	    */
	/************************************/

	for(i=1;i<=argc-1;i++){
		if(strcmp(argv[i],"-p")==0){   /* port number  */
			if(i+1<=argc-1){
				if(isint(argv[i+1])){
					port_number=atoi(argv[i+1]);	
					printf("port number is %i\n",port_number);
				}else{
					printf("bad argument\n");
					return 0;
				}
			}else{
				printf("bad argument\n");
				return 0;
			}
		}else if(strcmp(argv[i],"-d")==0){ /* debug flag */
			debug=1;
		}
	}

	/************************************/
	/*  print setting	            */
	/************************************/
	printf("WINDOW_TYPE\t=\t%s\n",num_to_str("WINDOW_TYPE",WINDOW_TYPE));
	printf("GRAPH_WINDOW\t=\t%s\n",one_to_yes(GRAPH_WINDOW));
	printf("RAND_TYPE\t=\t%i\n",RAND_TYPE);
	printf("RULE_KAKUMEI\t=\t%s\n",one_to_yes(RULE_KAKUMEI));
	printf("RULE_SHIBARI\t=\t%s\n",one_to_yes(RULE_SHIBARI));
	printf("RULE_KINSOKU\t=\t%s\n",one_to_yes(RULE_KINSOKU));
	printf("RULE_KAIDAN\t=\t%s\n",one_to_yes(RULE_KAIDAN));
	printf("RULE_CHANGE\t=\t%s\n",one_to_yes(RULE_CHANGE));
	printf("RULE_5TOBI\t=\t%s\n",one_to_yes(RULE_5TOBI));
	printf("RULE_6REVERS\t=\t%s\n",one_to_yes(RULE_6REVERS));
	printf("RULE_8GIRI\t=\t%s\n",one_to_yes(RULE_8GIRI));
	printf("RULE_11BACK\t=\t%s\n",one_to_yes(RULE_11BACK));
	printf("RULE_SEKIGAE\t=\t%s\n",one_to_yes(RULE_SEKIGAE));
	printf("RULE_SEKIGAE_NUM\t=\t%i\n",RULE_SEKIGAE_NUM);
	printf("GAME_NUMBER\t=\t%i\n",GAME_NUMBER);
	printf("FLASH_MIBUN_NUMBER\t=\t%i\n",FLASH_MIBUN_NUMBER);
	printf("GAME_PORT\t=\t%i\n",port_number);


	/************************************/
	/*  initialize random seed	    */
	/************************************/

	//srand((unsigned)time(NULL));
	tn_rand_init((unsigned long)time(NULL),RAND_TYPE);
	//tn_rand_init((unsigned long)3,RAND_TYPE);

	/************************************/
	/*  setting for X11   		    */
	/************************************/
	if(WINDOW_TYPE==0 || WINDOW_TYPE==1 || WINDOW_TYPE==3 || WINDOW_TYPE==4 || WINDOW_TYPE==5 || WINDOW_TYPE==6 ||GRAPH_WINDOW==1 || GRAPH_WINDOW==2){
		display = XOpenDisplay(NULL);
		root  = DefaultRootWindow( display );
		screen = DefaultScreen( display );
		background = WhitePixel(display, 0); //set windows
		foreground = BlackPixel(display, 0);
		attr.backing_store = Always; //backing store
	}

	if(WINDOW_TYPE==0 || WINDOW_TYPE==1 || WINDOW_TYPE==3 || WINDOW_TYPE==4 || WINDOW_TYPE==5 ||WINDOW_TYPE==6){
		switch(WINDOW_TYPE){
			case 0: // small
				win_info.card_width=38;
				win_info.card_height=57;
				win_info.card_offset_x=5+win_info.card_width*2;
				win_info.card_offset_y=30;
				win_info.window_width=586;
				win_info.window_height=467;
				win_info.mibun_width=win_info.card_width;
				win_info.mibun_height=win_info.card_height*8;
				win_info.control_width=50;
				win_info.control_height=200;

				win_info.on_stage_card_x_offset=5;
				win_info.on_stage_card_y_offset=20;
				win_info.on_stage_player_size=78;
				win_info.on_stage_player_x_offset=0;
				win_info.on_stage_player_y_offset=99;
				strcpy(XPM_CARD,"card_s.xpm");
				strcpy(XPM_TEFUDA,"tefuda_s.xpm");
				strcpy(XPM_TEFUDA2,"tefuda2_s.xpm");
				strcpy(XPM_MIBUN,"mibun_s.xpm");
				break;
			case 1: // big
				win_info.card_width=77;
				win_info.card_height=114;
				win_info.card_offset_x=5+win_info.card_width*2;
				win_info.card_offset_y=30;
				win_info.window_width=1192;
				win_info.window_height=933;
				win_info.mibun_width=win_info.card_width;
				win_info.mibun_height=win_info.card_height*8;
				win_info.control_width=50;
				win_info.control_height=200;

				win_info.on_stage_card_x_offset=5;
				win_info.on_stage_card_y_offset=41;
				win_info.on_stage_player_size=155;
				win_info.on_stage_player_x_offset=0;
				win_info.on_stage_player_y_offset=198;
				strcpy(XPM_CARD,"card_b.xpm");
				strcpy(XPM_TEFUDA,"tefuda_b.xpm");
				strcpy(XPM_TEFUDA2,"tefuda2_b.xpm");
				strcpy(XPM_MIBUN,"mibun_b.xpm");
				break;
			case 4: // middle
				win_info.card_width=45;
				win_info.card_height=67;
				win_info.card_offset_x=5+win_info.card_width*2;
				win_info.card_offset_y=30;
				win_info.window_width=703;
				win_info.window_height=550;
				win_info.mibun_width=win_info.card_width;
				win_info.mibun_height=win_info.card_height*8;
				win_info.control_width=50;
				win_info.control_height=200;

				win_info.on_stage_card_x_offset=5;
				win_info.on_stage_card_y_offset=24;
				win_info.on_stage_player_size=24+win_info.card_height;
				win_info.on_stage_player_x_offset=0;
				win_info.on_stage_player_y_offset=136;
				strcpy(XPM_CARD,"card_m.xpm");
				strcpy(XPM_TEFUDA,"tefuda_m.xpm");
				strcpy(XPM_TEFUDA2,"tefuda2_m.xpm");
				strcpy(XPM_MIBUN,"mibun_m.xpm");
				break;
			case 3: // middle extend
				win_info.card_width=45;
				win_info.card_height=67;
				win_info.card_offset_x=5+win_info.card_width*2;
				win_info.card_offset_y=30;
				win_info.window_width=703;
				win_info.window_height=570;
				win_info.mibun_width=win_info.card_width;
				win_info.mibun_height=win_info.card_height*8;
				win_info.control_width=50;
				win_info.control_height=200;

				win_info.on_stage_card_x_offset=5;
				win_info.on_stage_card_y_offset=44;
				win_info.on_stage_player_size=24+win_info.card_height;
				win_info.on_stage_player_x_offset=0;
				win_info.on_stage_player_y_offset=136;
				strcpy(XPM_CARD,"card_me.xpm");
				strcpy(XPM_TEFUDA,"tefuda_me.xpm");
				strcpy(XPM_TEFUDA2,"tefuda2_me.xpm");
				strcpy(XPM_MIBUN,"mibun_me.xpm");
				break;
			case 5: // middle xga
				win_info.card_width=56;
				win_info.card_height=84;
				win_info.card_offset_x=5+win_info.card_width*2;
				win_info.card_offset_y=30;
				win_info.window_width=950;
				win_info.window_height=720;
				win_info.mibun_width=win_info.card_width;
				win_info.mibun_height=win_info.card_height*8;
				win_info.control_width=50;
				win_info.control_height=200;

				win_info.on_stage_card_x_offset=5;
				win_info.on_stage_card_y_offset=52;
				win_info.on_stage_player_size=32+win_info.card_height;
				win_info.on_stage_player_x_offset=100;
				win_info.on_stage_player_y_offset=168;
				sprintf(XPM_CARD,"%scard_l.xpm",XPM_PATH);
				sprintf(XPM_TEFUDA,"%stefuda_l.xpm",XPM_PATH);
				sprintf(XPM_TEFUDA2,"%stefuda2_l.xpm",XPM_PATH);
				sprintf(XPM_MIBUN,"%smibun_l.xpm",XPM_PATH);

				//strcpy(XPM_CARD,"card_l.xpm");
				//strcpy(XPM_TEFUDA,"tefuda_l.xpm");
				//strcpy(XPM_TEFUDA2,"tefuda2_l.xpm");
				//strcpy(XPM_MIBUN,"mibun_l.xpm");
				break;
			case 6: // NEW_NORMAL
				win_info.window_width=900;
				win_info.window_height=700;
				win_info.control_width=50;
				win_info.control_height=250;
				win_info.card_width=70;
				win_info.card_height=80;
				win_info.mibun_width=150;
				win_info.mibun_height=60;
				win_info.card_offset_x=5+win_info.card_width*2;
				win_info.card_offset_y=30;
				win_info.on_stage_card_x_offset=5;
				win_info.on_stage_card_y_offset=52;
				win_info.on_stage_player_size=32+win_info.card_height;
				win_info.on_stage_player_x_offset=100;
				win_info.on_stage_player_y_offset=168;
				sprintf(XPM_CARD,"%scard_newn.xpm",XPM_PATH);
				sprintf(XPM_MIBUN,"%smibun_newn.xpm",XPM_PATH);
				sprintf(XPM_TEFUDA,"%stefuda_newn.xpm",XPM_PATH);
				sprintf(XPM_TEFUDA2,"%stefuda2_newn.xpm",XPM_PATH);
				break;
		}	
		
		win_info.root=root;
		win_info.screen=screen;
		(win_info).display=display;
	
		win_info.window_type=WINDOW_TYPE;
		
		win_info.screen0 = XCreatePixmap(win_info.display, // tmp screen
			DefaultRootWindow(win_info.display),
			win_info.window_width, win_info.window_height,
			DefaultDepth(win_info.display,win_info.screen ));
		win_info.screen1 = XCreatePixmap(win_info.display, // screen1 : background bitmap in the case of normal
			DefaultRootWindow(win_info.display),
			win_info.window_width, win_info.window_height,
			DefaultDepth(win_info.display, win_info.screen ));
		win_info.screen2 = XCreatePixmap(win_info.display, // screen2 : card bitmap
			DefaultRootWindow(win_info.display),
			win_info.window_width, win_info.window_height,
			DefaultDepth(win_info.display, win_info.screen));
		win_info.screen3 = XCreateSimpleWindow(win_info.display, // display screen
			DefaultRootWindow(win_info.display),
			0, 0, win_info.window_width, win_info.window_height,
			0, 0, win_info.background);
		XSetStandardProperties(win_info.display, win_info.screen3, 
			"daihinmin", "daihugou",
			None, argv, argc, NULL);
		win_info.screen4 = XCreatePixmap(win_info.display, // screen4 : background bitmap in the case of kakumei
			DefaultRootWindow(win_info.display),
			win_info.window_width, win_info.window_height,
			DefaultDepth(win_info.display, win_info.screen ));
		if(WINDOW_TYPE==6){
			win_info.screen5 = XCreatePixmap(win_info.display, // tmp screen
				DefaultRootWindow(win_info.display),
				1201,960,
				DefaultDepth(win_info.display,win_info.screen ));
			win_info.screen6 = XCreatePixmap(win_info.display, // tmp screen
				DefaultRootWindow(win_info.display),
				1500,800,
				DefaultDepth(win_info.display,win_info.screen ));
		}
		win_info.control_win = XCreateSimpleWindow(win_info.display, // control_win : control window bitmsp
			DefaultRootWindow(win_info.display),
			0, 0, win_info.control_width, win_info.control_height,
			0, 0, win_info.background);
		XSetStandardProperties(win_info.display, win_info.control_win, 
			"control", "daihugou",
			None, argv, argc, NULL);
		win_info.mibun_win = XCreatePixmap(win_info.display, // mibun bitmap
			DefaultRootWindow(win_info.display),
			win_info.mibun_width, win_info.mibun_height,
			DefaultDepth(win_info.display, win_info.screen ));

		XSelectInput(win_info.display, win_info.control_win, ExposureMask | ButtonPressMask);

		XChangeWindowAttributes( win_info.display, win_info.screen3, CWBackingStore, &attr );
		XChangeWindowAttributes( win_info.display, win_info.control_win, CWBackingStore, &attr );

		XMapWindow(win_info.display, win_info.screen3);
		XMapWindow(win_info.display, win_info.control_win);
		
		win_info.background = WhitePixel(display, 0); //set windows
		win_info.foreground = BlackPixel(display, 0);

		win_info.gc = XCreateGC(win_info.display, win_info.screen0, 0, 0);
		win_info.font=XLoadFont(win_info.display,"-adobe-times-bold-i-normal--*-240-*-*-p-*-iso8859-1");
		XSetFont(win_info.display, win_info.gc, win_info.font);
		XSetForeground(win_info.display,win_info.gc,win_info.foreground);
		XSetBackground(win_info.display,win_info.gc,win_info.background);

		if(WINDOW_TYPE==0 || WINDOW_TYPE==1 || WINDOW_TYPE==3 || WINDOW_TYPE==4 || WINDOW_TYPE==5 ){
			err = XpmReadFileToPixmap(win_info.display, win_info.screen2,XPM_CARD,&win_info.pix,NULL,NULL);
			XCopyArea(win_info.display, win_info.pix, win_info.screen2, win_info.gc, 0,0,win_info.window_width,win_info.window_height,0,0); // screen2 : card bitmap
			err = XpmReadFileToPixmap(win_info.display,win_info.screen1,XPM_TEFUDA,&win_info.pix,NULL,NULL);
			XCopyArea(win_info.display, win_info.pix, win_info.screen1, win_info.gc, 0,0,win_info.window_width,win_info.window_height,0,0); // screen1 : background bitmap in the case of normal
			XCopyArea(win_info.display, win_info.pix, win_info.screen0, win_info.gc, 0,0,win_info.window_width,win_info.window_height,0,0); // screen0 : tmp screen
			err = XpmReadFileToPixmap(win_info.display,win_info.screen1,XPM_TEFUDA2,&win_info.pix,NULL,NULL);  
			XCopyArea(win_info.display, win_info.pix, win_info.screen4, win_info.gc, 0,0,win_info.window_width,win_info.window_height,0,0); // screen4 : background bitmap in the case of kakumei
			err = XpmReadFileToPixmap(win_info.display,win_info.mibun_win,XPM_MIBUN,&win_info.pix,NULL,NULL);
			XCopyArea(win_info.display, win_info.pix, win_info.mibun_win, win_info.gc, 0,0,win_info.mibun_width,win_info.mibun_height,0,0); // mibun   : mibun bitmap
		}
		if(WINDOW_TYPE==6){
			err = XpmReadFileToPixmap(win_info.display,win_info.screen1,XPM_TEFUDA,&win_info.pix,NULL,NULL);
			XCopyArea(win_info.display, win_info.pix, win_info.screen1, win_info.gc, 0,0,win_info.window_width,win_info.window_height,0,0); // screen1 : background bitmap in the case of normal
			XCopyArea(win_info.display, win_info.pix, win_info.screen0, win_info.gc, 0,0,win_info.window_width,win_info.window_height,0,0); // screen0 : tmp screen
			err = XpmReadFileToPixmap(win_info.display,win_info.screen1,XPM_TEFUDA2,&win_info.pix,NULL,NULL);  
			XCopyArea(win_info.display, win_info.pix, win_info.screen4, win_info.gc, 0,0,win_info.window_width,win_info.window_height,0,0); // screen4 : background bitmap in the case of kakumei

			err = XpmReadFileToPixmap(win_info.display, win_info.screen5,XPM_CARD,&win_info.pix,NULL,NULL);
			XCopyArea(win_info.display, win_info.pix, win_info.screen5, win_info.gc, 0,0,1201,960,0,0); // screen5 : card bitmap
			err = XpmReadFileToPixmap(win_info.display, win_info.screen6,XPM_MIBUN,&win_info.pix,NULL,NULL);
			XCopyArea(win_info.display, win_info.pix, win_info.screen6, win_info.gc, 0,0,1500,800,0,0); // screen6 : win
		}

		XCopyArea(win_info.display, win_info.screen0, win_info.screen3, win_info.gc, 0,0,win_info.window_width, win_info.window_height, 0,0); // delete old card on stage

                XSetForeground(win_info.display,win_info.gc,win_info.background); // clean pixmap
                XSetBackground(win_info.display,win_info.gc,win_info.background);
                XFillRectangle(win_info.display, win_info.control_win, win_info.gc, 0,0, win_info.control_width, win_info.control_height);
                XSetForeground(win_info.display,win_info.gc,win_info.foreground); // clean pixmap
		for(i=1;i<=3;i++){
			XDrawLine(win_info.display, win_info.control_win, win_info.gc, 
				0, i*50, win_info.control_width, i*50);
		}

		XFlush( display );

	} // fi (WINDOW_TYPE == 0 or 1 or 2 or 3 or 4 or 5)

	if(GRAPH_WINDOW==1 || GRAPH_WINDOW==2){

		switch(GRAPH_WINDOW){
			case 1:
				g_win_info.window_width=1100;
				g_win_info.window_height=410;
				g_win_info.x_scale=20;
				g_win_info.y_scale=4;
				g_win_info.x_border=g_win_info.x_scale*50;
				g_win_info.y_border=g_win_info.y_scale*100;
				break;
			case 2:
				g_win_info.window_width=900;
				g_win_info.window_height=700;
				g_win_info.x_scale=20;
				g_win_info.y_scale=6;
				g_win_info.x_border=g_win_info.x_scale*50/4*3;
				//g_win_info.y_border=g_win_info.y_scale*100;
				g_win_info.y_border=600;
				break;
		}
		g_win_info.y_scale2=1.5/2/2/2;

		g_win_info.background = WhitePixel(display, 0); //set windows
		g_win_info.foreground = BlackPixel(display, 0);
		g_win_info.root=root;
		g_win_info.screen=screen;
		(g_win_info).display=display;

		g_win_info.window = XCreateSimpleWindow(g_win_info.display,
			DefaultRootWindow(g_win_info.display),
			0, 0, g_win_info.window_width, g_win_info.window_height,
			0, 0, g_win_info.background);
		g_win_info.pixmap = XCreatePixmap(g_win_info.display, 
			DefaultRootWindow(g_win_info.display),
			g_win_info.window_width, g_win_info.window_height,
			DefaultDepth(g_win_info.display,g_win_info.screen));
		g_win_info.pixmap2 = XCreatePixmap(g_win_info.display, 
			DefaultRootWindow(g_win_info.display),
			g_win_info.window_width, g_win_info.window_height,
			DefaultDepth(g_win_info.display,g_win_info.screen));
		XSetStandardProperties(g_win_info.display, g_win_info.window, 
			"control", "graph",
			None, argv, argc, NULL);
		XChangeWindowAttributes( g_win_info.display, g_win_info.window, CWBackingStore, &attr );

		gc = XCreateGC(display, g_win_info.window, 0, 0);
		g_win_info.gc=gc;
		g_win_info.font=XLoadFont(g_win_info.display,"-adobe-times-bold-i-normal--*-240-*-*-p-*-iso8859-1");
		//g_win_info.font=XLoadFont(g_win_info.display,"9x15bold");
		XSetFont(g_win_info.display, g_win_info.gc, g_win_info.font);


		XMapWindow(g_win_info.display,g_win_info.window);
		graph_initialize(&g_win_info, player_name);
		graph_initialize2(&g_win_info, player_name);
		XFlush( g_win_info.display );


	}//fi GRAPH_WINDOW

	/************************************/
	/*  file open for log	       */
	/************************************/

	fp=fopen("debug.dat","w");
	fp2=fopen("debug2.dat","w");

	/********************************/
	/* setting for client/server	*/
	/*  make soket for each client  */
	/********************************/
	if ((sockfd = socket(PF_INET, SOCK_STREAM, 0)) < 0) {
		perror("client: socket");
		exit(1);
	}
	bzero((char *) &wait_addr, sizeof(wait_addr));
	wait_addr.sin_family = PF_INET;
	wait_addr.sin_addr.s_addr = htons(INADDR_ANY);
	wait_addr.sin_port = htons(port_number);

	i = 1;
	j = sizeof(i);
	if (setsockopt(sockfd, SOL_SOCKET, SO_REUSEADDR, (char *)&i, j) < 0){
		perror("setsockopt");
	}

	if (bind(sockfd,(struct sockaddr *)&wait_addr,sizeof(wait_addr))<0) {
		perror("reader: bind");
		exit(1);
	}
	if (listen(sockfd, 1) < 0) {
		perror("reader: listen");
		close(sockfd);
		exit(1);
	}
	for(i=0;i<=4;i++){
		printf("now waiting %i \n", i);
		client_len[i]=sizeof(client_addr[i]);
		if((client_sockfd[i]=accept(sockfd,(struct sockaddr *)&client_addr[i],&client_len[i])) < 0 ){
			perror("reader: accept");
			exit(1);
		};
		FD_ZERO(&org_target_fds);
		FD_SET(client_sockfd[i], &org_target_fds);
		memcpy(&target_fds, &org_target_fds, sizeof(org_target_fds));
		waitval.tv_sec  = 2;
		waitval.tv_usec = 500;
		switch(select(50,&target_fds,NULL,NULL,&waitval)){
			case -1:
				printf("protocol_version: NONE\n");
				exit(1);
			case 0: /* time out */
				protocol_version=20060;
				printf("protocol_version: 2006a\n");
				break;
			default: /* connect from client */
				tn_card_read(client_sockfd[i], work_card , protocol_version);
				protocol_version=work_card[0][0];
				printf("protocol_version: %d\n",work_card[0][0]);
				for(j=0;j<=8;j++){
					player_name[i][j]=(char)work_card[1][j];
				}
				player_name[i][9]='\0';
				printf("NAME: %s\n", player_name[i]);
				break;
		}
		tn_int_write(client_sockfd[i], i , protocol_version);
		printf("accepted from %s \n",inet_ntoa(client_addr[i].sin_addr));
	} 
	/* end of socket setting */


	/*********************/
	/* initialize values */
	/*********************/

	for(i=0;i<=4;i++){
		mibun[i]=i;	
	}
	for(i=0;i<=4;i++){
		sekijun[i]=i;	
	}

	// player statistics
	
	for(i=0;i<=4;i++){
		bzero(&ps[i], sizeof(struct playerStatistics));
	}

	/**************/
	/* game start */
	/**************/
	for(now_number_of_games=1;now_number_of_games<=GAME_NUMBER;now_number_of_games++){
		// shuffle all cards for 5 players
		bzero(stage_card,sizeof(stage_card));
		bzero(old_stage_card,sizeof(old_stage_card));
		bzero(players_card,sizeof(players_card));
		bzero(goal_flag,sizeof(goal_flag));

		if(RULE_SEKIGAE!=0){ // decide sekigae
			tn_sekigae(now_number_of_games,sekijun,RULE_SEKIGAE,RULE_SEKIGAE_NUM, RAND_TYPE);
			printf("sekigae done\n");
			for(i=0;i<=4;i++){
				printf("-> %i\n",sekijun[i]);
			}
		}

		if(((now_number_of_games-1 )% FLASH_MIBUN_NUMBER)==0){ // decide first member for shuffle
			//shuffle_card(rand()%5,players_card);
			shuffle_card((int)(tn_rand_gen(RAND_TYPE)*5),players_card, RAND_TYPE, sekijun);
		}else{
			shuffle_card(mibun[0],players_card , RAND_TYPE, sekijun);
		}

		for(i=0;i<=4;i++){
			work_card[6][mibun[i]+5]=i;
		}
		for(i=0;i<=4;i++){ // initialize each table[5]
			bcopy(work_card[5], players_card[i][5], 2*sizeof(work_card[5]));
		}

		if(WINDOW_TYPE==0 || WINDOW_TYPE==1 || WINDOW_TYPE==3 || WINDOW_TYPE==4){
			for(i=0;i<=4;i++){ //  all user is heimin
				if(((now_number_of_games-1 )% FLASH_MIBUN_NUMBER)==0){
					XCopyArea(win_info.display,win_info.mibun_win,win_info.screen1,win_info.gc,0,win_info.card_height*3,win_info.card_width,win_info.card_height,win_info.on_stage_card_x_offset,win_info.on_stage_player_size*i+win_info.on_stage_player_y_offset);	
					XCopyArea(win_info.display,win_info.mibun_win,win_info.screen4,win_info.gc,0,win_info.card_height*3,win_info.card_width,win_info.card_height,win_info.on_stage_card_x_offset,win_info.on_stage_player_size*i+win_info.on_stage_player_y_offset);	
				}else{
					XCopyArea(win_info.display,win_info.mibun_win,win_info.screen1,win_info.gc,0,win_info.card_height*(players_card[0][6][5+i]+1),win_info.card_width,win_info.card_height,win_info.on_stage_card_x_offset,win_info.on_stage_player_size*i+win_info.on_stage_player_y_offset);	
					XCopyArea(win_info.display,win_info.mibun_win,win_info.screen4,win_info.gc,0,win_info.card_height*(players_card[0][6][5+i]+1),win_info.card_width,win_info.card_height,win_info.on_stage_card_x_offset,win_info.on_stage_player_size*i+win_info.on_stage_player_y_offset);	
				}
			}
		}
		if(WINDOW_TYPE==5){
			for(i=0;i<=4;i++){ //  all user is heimin
				if(((now_number_of_games-1 )% FLASH_MIBUN_NUMBER)==0){
					XCopyArea(win_info.display,win_info.mibun_win,win_info.screen1,win_info.gc,0,win_info.card_height*3,win_info.card_width,win_info.card_height,win_info.on_stage_card_x_offset+win_info.on_stage_player_x_offset,win_info.on_stage_player_size*i+win_info.on_stage_player_y_offset);	
					XCopyArea(win_info.display,win_info.mibun_win,win_info.screen4,win_info.gc,0,win_info.card_height*3,win_info.card_width,win_info.card_height,win_info.on_stage_card_x_offset+win_info.on_stage_player_x_offset,win_info.on_stage_player_size*i+win_info.on_stage_player_y_offset);	
				}else{
					XCopyArea(win_info.display,win_info.mibun_win,win_info.screen1,win_info.gc,0,win_info.card_height*(players_card[0][6][5+i]+1),win_info.card_width,win_info.card_height,win_info.on_stage_card_x_offset+win_info.on_stage_player_x_offset,win_info.on_stage_player_size*get_seat(sekijun,i)+win_info.on_stage_player_y_offset);	
					XCopyArea(win_info.display,win_info.mibun_win,win_info.screen4,win_info.gc,0,win_info.card_height*(players_card[0][6][5+i]+1),win_info.card_width,win_info.card_height,win_info.on_stage_card_x_offset+win_info.on_stage_player_x_offset,win_info.on_stage_player_size*get_seat(sekijun,i)+win_info.on_stage_player_y_offset);	
				}
			}
		}
		if(WINDOW_TYPE==6){
				if(((now_number_of_games-1 )% FLASH_MIBUN_NUMBER)==0){
					XCopyArea(win_info.display,win_info.screen5,win_info.screen1,win_info.gc,300,0,150,60,130,5);	
					XCopyArea(win_info.display,win_info.screen5,win_info.screen1,win_info.gc,300,0,150,60,430,5);	
					XCopyArea(win_info.display,win_info.screen5,win_info.screen1,win_info.gc,300,0,150,60,730,5);	
					XCopyArea(win_info.display,win_info.screen5,win_info.screen1,win_info.gc,300,0,150,60,130,355);	
					XCopyArea(win_info.display,win_info.screen5,win_info.screen1,win_info.gc,300,0,150,60,730,355);	
					XCopyArea(win_info.display,win_info.screen5,win_info.screen4,win_info.gc,300,0,150,60,130,5);	
					XCopyArea(win_info.display,win_info.screen5,win_info.screen4,win_info.gc,300,0,150,60,430,5);	
					XCopyArea(win_info.display,win_info.screen5,win_info.screen4,win_info.gc,300,0,150,60,730,5);	
					XCopyArea(win_info.display,win_info.screen5,win_info.screen4,win_info.gc,300,0,150,60,130,355);	
					XCopyArea(win_info.display,win_info.screen5,win_info.screen4,win_info.gc,300,0,150,60,730,355);	

				}else{
					XCopyArea(win_info.display,win_info.screen5,win_info.screen1,win_info.gc,(4-get_seat(mibun,sekijun[0]))*150,0,150,60,130,5);	
					XCopyArea(win_info.display,win_info.screen5,win_info.screen1,win_info.gc,(4-get_seat(mibun,sekijun[1]))*150,0,150,60,430,5);	
					XCopyArea(win_info.display,win_info.screen5,win_info.screen1,win_info.gc,(4-get_seat(mibun,sekijun[2]))*150,0,150,60,730,5);	
					XCopyArea(win_info.display,win_info.screen5,win_info.screen1,win_info.gc,(4-get_seat(mibun,sekijun[3]))*150,0,150,60,130,355);	
					XCopyArea(win_info.display,win_info.screen5,win_info.screen1,win_info.gc,(4-get_seat(mibun,sekijun[4]))*150,0,150,60,730,355);	
					XCopyArea(win_info.display,win_info.screen5,win_info.screen4,win_info.gc,(4-get_seat(mibun,sekijun[0]))*150,0,150,60,130,5);	
					XCopyArea(win_info.display,win_info.screen5,win_info.screen4,win_info.gc,(4-get_seat(mibun,sekijun[1]))*150,0,150,60,430,5);	
					XCopyArea(win_info.display,win_info.screen5,win_info.screen4,win_info.gc,(4-get_seat(mibun,sekijun[2]))*150,0,150,60,730,5);	
					XCopyArea(win_info.display,win_info.screen5,win_info.screen4,win_info.gc,(4-get_seat(mibun,sekijun[3]))*150,0,150,60,130,355);	
					XCopyArea(win_info.display,win_info.screen5,win_info.screen4,win_info.gc,(4-get_seat(mibun,sekijun[4]))*150,0,150,60,730,355);	
				}
		}
	
		if(WINDOW_TYPE==0 || WINDOW_TYPE==1 || WINDOW_TYPE==3 || WINDOW_TYPE==4 || WINDOW_TYPE==5 || WINDOW_TYPE==6){
			if(flag_wait_type!=3){
				if(WINDOW_TYPE==6){
					tmp=stage_card[5][0];
					if(((now_number_of_games % FLASH_MIBUN_NUMBER)!=1)&&(RULE_CHANGE==1)){ 
						stage_card[5][0]=1;
					}
					initialize_window2(stage_card, old_stage_card, players_card,&win_info,now_pass, player_name, sekijun, &accept_flag, mibun);
					stage_card[5][0]=tmp;
				}else{
					initialize_window(stage_card, old_stage_card, players_card,&win_info,now_pass, player_name, sekijun);
				}
				if((now_number_of_games % FLASH_MIBUN_NUMBER)!=1){
					if(WINDOW_TYPE==5){
						XCopyArea(win_info.display,win_info.screen1,win_info.screen3,win_info.gc,
							win_info.on_stage_card_x_offset+4*win_info.card_width,win_info.on_stage_card_y_offset-20, 
							4*win_info.card_width,win_info.card_height,
							win_info.on_stage_card_x_offset+4*win_info.card_width,win_info.on_stage_card_y_offset-20);	
						XCopyArea(win_info.display,win_info.screen2,win_info.screen3,win_info.gc,
							win_info.card_width*2,win_info.card_height*4,3*win_info.card_width,win_info.card_height,
							win_info.on_stage_card_x_offset+4*win_info.card_width,win_info.on_stage_card_y_offset-20);	
					}
				}
			}
			wait_control(&win_info,&g_win_info, &flag_wait_type, &graph_mode, &g_flag);
		}
		if(debug){printf("shuffle OK\n");} // DEBUG

		/************************************/
		/* distribute cards to each players */
		/************************************/
		//convert status to table
		work_card[5][0] = 0;  
		work_card[5][1] = 0;
		work_card[5][2] = 0;
		work_card[5][3] = 6; 
		work_card[5][4] = 0;
		work_card[5][5] = 0; 
		work_card[5][6] = 0;
		work_card[5][7] = 0;
		work_card[6][0] = count_card_num_r(players_card[0], &error);
		work_card[6][1] = count_card_num_r(players_card[1], &error);
		work_card[6][2] = count_card_num_r(players_card[2], &error);
		work_card[6][3] = count_card_num_r(players_card[3], &error);
		work_card[6][4] = count_card_num_r(players_card[4], &error);
		if(((now_number_of_games-1 )% FLASH_MIBUN_NUMBER)==0){
			for(i=0;i<=4;i++){
				work_card[6][i+5]=0;
			}
		}else{
			for(i=0;i<=4;i++){
				work_card[6][mibun[i]+5]=i;
			}
		}
		for(i=10;i<=14;i++){
			work_card[6][i]=sekijun[i-10];
		}
		for(i=0;i<=4;i++){
			bcopy(work_card[5], players_card[i][5], 2*sizeof(work_card[5]));
		}
		// distribute

	
		// search strong card in hinmin and data setting
		if(((now_number_of_games % FLASH_MIBUN_NUMBER)!=1)&&(RULE_CHANGE==1)){ 
			for(i=0;i<=4;i++){  // data setting for all player
				players_card[mibun[i]][5][0]=1;
				players_card[mibun[i]][5][1]=2-i;
			}
			bcopy(players_card[mibun[3]], tmp_card[0], sizeof(players_card[mibun[3]]));
			bcopy(players_card[mibun[4]], tmp_card[1], sizeof(players_card[mibun[4]]));
			for(k=2;k>=1;k--){ //search strong card
				trans_strong_card(players_card[mibun[2+k]],players_card[mibun[2-k]],k);
			}
		}else{  // no change because of it's a 1st game
			for(i=0;i<=4;i++){  // data setting for all player
				players_card[mibun[i]][5][0]=1;  
				players_card[mibun[i]][5][1]=0;
			}
			bcopy(players_card[mibun[3]], tmp_card[0], sizeof(players_card[mibun[3]]));
			bcopy(players_card[mibun[4]], tmp_card[1], sizeof(players_card[mibun[4]]));

		}
		if(debug){printf("hinmin and daihinmin card change is done\n");} //DEBUG
		tn_card_write(client_sockfd[mibun[0]],players_card[mibun[0]],protocol_version); // tuuchi
		tn_card_write(client_sockfd[mibun[1]],players_card[mibun[1]],protocol_version);
		tn_card_write(client_sockfd[mibun[2]],players_card[mibun[2]],protocol_version);
		//tn_card_write(client_sockfd[mibun[3]],players_card[mibun[3]],protocol_version);
		//tn_card_write(client_sockfd[mibun[4]],players_card[mibun[4]],protocol_version);
		tn_card_write(client_sockfd[mibun[3]],tmp_card[0],protocol_version);
		tn_card_write(client_sockfd[mibun[4]],tmp_card[1],protocol_version);

				// change faise of daihugou start
		if(((now_number_of_games % FLASH_MIBUN_NUMBER)!=1)&&(RULE_CHANGE==1)){ // change card
			tn_card_read(client_sockfd[mibun[0]],work_card,protocol_version); // uketori
			error=0;
			error=count_card_num(work_card, &number_of_card); // number check
			if((check_include_card(players_card[mibun[0]],work_card)==0)&&(number_of_card==2)&&(error==0)){
				if(debug){printf("change card - OK \n");}
				trans_work_card(players_card[mibun[0]],players_card[mibun[4]],work_card);
			}else{
				if(debug){printf("change card - fault \n");}
				trans_strong_card(players_card[mibun[0]],players_card[mibun[4]],2);
			} // fi
		}// fi
		if(debug){printf("change daihugou - OK\n");} //DEBUG

				// change faise of hugou start
		if(((now_number_of_games % FLASH_MIBUN_NUMBER)!=1)&&(RULE_CHANGE==1)){ // change card
			tn_card_read(client_sockfd[mibun[1]],work_card,protocol_version);
			error=0;
			error=count_card_num(work_card, &number_of_card); // number check
			if((check_include_card(players_card[mibun[1]],work_card)==0)&&(number_of_card==1)){
				if(debug){printf("change card - OK \n");}
				trans_work_card(players_card[mibun[1]],players_card[mibun[3]],work_card);
			}else{
				if(debug){printf("change card - fault \n");}
				trans_strong_card(players_card[mibun[1]],players_card[mibun[3]],1);
			} // fi
		}// fi
		if(debug){printf("change hugou - OK\n");} //DEBUG
		if(debug){printf("end of distribute\n");} //DEBUG
		/* end of distribute */


		/****************************/
		/* each players phase start */
		/****************************/

		status_of_submitted_card[0]=-1;    // flash of stage
		status_of_stages_card[0]=-1;    // flash of stage
		status_of_stages_card[1]=-1;   
		status_of_stages_card[2]=-1;   
		status_of_stages_card[3]=-1;   
		bzero(stage_card,sizeof(stage_card));
		bzero(work_card,sizeof(work_card));
		now_number_of_goal=0;
		now_player=search_card(players_card,2,1);    // search a "dia 3"
		last_player=now_player; 

		now_muki=1;
		now_kakumei=0;
		now_11back=0;
		now_jocker_single=0;   // flag of jocker
		now_shibari=0;   // flag of jocker
		number_renzoku_pass=0;
	
		flash=0;   // flash
		bzero(goal_flag,sizeof(goal_flag));      // flags of person  
		bzero(now_pass,sizeof(now_pass));       // flags of passed person

		// player statistics
		count_turn = 0;

		while(now_number_of_goal<=3){

			//convert status to table
			work_card[5][0] = 0;  
			work_card[5][1] = 0;
			work_card[5][2] = 0;
			work_card[5][3] = now_player; 
			work_card[5][4] = (status_of_stages_card[0]==-1);
			work_card[5][5] = now_11back; 
			work_card[5][6] = now_kakumei;
			work_card[5][7] = now_shibari;
			work_card[6][0] = count_card_num_r(players_card[0], &error);
			work_card[6][1] = count_card_num_r(players_card[1], &error);
			work_card[6][2] = count_card_num_r(players_card[2], &error);
			work_card[6][3] = count_card_num_r(players_card[3], &error);
			work_card[6][4] = count_card_num_r(players_card[4], &error);
			if(((now_number_of_games-1 )% FLASH_MIBUN_NUMBER)==0){
				for(i=0;i<=4;i++){
					work_card[6][i+5]=0;
				}
			}else{
				for(i=0;i<=4;i++){
					work_card[6][mibun[i]+5]=i;
				}
			}
			for(i=10;i<=14;i++){
				work_card[6][i]=sekijun[i-10];
			}
			for(i=0;i<=4;i++){
				bcopy(work_card[5], players_card[i][5], 2*sizeof(work_card[5]));
			}
			bcopy(work_card[5], stage_card[5], 2*sizeof(work_card[5]));
			players_card[now_player][5][2]=1;

			if(WINDOW_TYPE==0 || WINDOW_TYPE==1 || WINDOW_TYPE==3 || WINDOW_TYPE==4 || WINDOW_TYPE==5 || WINDOW_TYPE==6){
				if(flag_wait_type!=3){
					if(WINDOW_TYPE==6){
						initialize_window2(stage_card, old_stage_card, players_card,&win_info,now_pass , player_name, sekijun, &accept_flag, mibun);
					}else{
						initialize_window(stage_card, old_stage_card, players_card,&win_info,now_pass , player_name, sekijun);
					}
				}
				wait_control(&win_info,&g_win_info, &flag_wait_type, &graph_mode, &g_flag);
			}

			if(debug){printf("To send prepared datas for each player now\n");} //DEBUG
			// prepare datas for each player
				// atode
			tn_card_write(client_sockfd[0],players_card[0],protocol_version);
			tn_card_write(client_sockfd[1],players_card[1],protocol_version);
			tn_card_write(client_sockfd[2],players_card[2],protocol_version);
			tn_card_write(client_sockfd[3],players_card[3],protocol_version);
			tn_card_write(client_sockfd[4],players_card[4],protocol_version);
			if(debug){printf("To send prepared date is done. \n");} //DEBUG

			if(debug){printf("Read from player %i\n", now_player);} //DEBUG
			tn_card_read(client_sockfd[now_player],work_card,protocol_version);
			if(debug){printf("accepted card is \n");} //DEBUG
			if(debug){print_player_card(work_card);} //DEBUG
			if(debug){printf("player card is \n");} //DEBUG
			if(debug){print_player_card(players_card[now_player]);} //DEBUG

			error=0; // error initialize

			error|=check_include_card(players_card[now_player],work_card);  // include check
			if(debug){printf("error01 = %i: include check\n",error);} //DEBUG

			error|=analyze_card(work_card,status_of_submitted_card,0); // analyze check.  the number of joker > 2 then error and etc
			if(debug){printf("error02 = %i: analyze check \n",error);} //DEBUG

			error|=compare_card_status(status_of_stages_card, status_of_submitted_card, (now_kakumei + now_11back)%2); // compare check
			if(debug){printf("error03 = %i: compare check \n",error);} //DEBUG
			if(debug){
				printf(
					"stage value    %i - number %i - suit %i - type %i \n",
					status_of_stages_card[0],status_of_stages_card[1],status_of_stages_card[2],status_of_stages_card[3]
				); //DEBUG
				printf(
					"submitte value %i - number %i - suit %i - type %i \n",
					status_of_submitted_card[0],status_of_submitted_card[1],status_of_submitted_card[2],status_of_submitted_card[3]
				); //DEBUG
			}

			count_turn++;
			if ((error==0)&&(status_of_stages_card[0]==-1)) {
				ps[now_player].getStage++;
				ps[now_player].cardStrength+=status_of_submitted_card[0];
				switch(status_of_submitted_card[3]) {
				case 1:
					break;
				case 2:
					ps[now_player].fukusuu++;
					break;
				case 3:
					ps[now_player].kaidan++;
					break;
				case 4:
					ps[now_player].jokerCnt++;
					ps[now_player].jokerTurnSum+=count_turn;
					break;
				default:
					break;
				}
			} else if ((error==0)) {
				switch(status_of_submitted_card[3]) {
				case 1:
					break;
				case 2:
					break;
				case 3:
					break;
				case 4:
					ps[now_player].jokerCnt++;
					ps[now_player].jokerTurnSum+=count_turn;
					break;
				default:
					break;
				}
			}

			if(RULE_SHIBARI){
				if((error==0)&&(status_of_submitted_card[3]!=4)){
					if(now_shibari==1){
						if(status_of_submitted_card[2]==status_of_stages_card[2]){

						}else{
							error|=1;
							if(debug){printf("error06 = %i: shibari \n", error);} //DEBUG
						}
					}else{
						ps[now_player].shibariCnt++;
						if(status_of_submitted_card[2]==status_of_stages_card[2]){
							now_shibari++;
							ps[now_player].shibari++;
						}else{
							now_shibari=0;
						}
					}
				}
			 } 

			if(now_jocker_single){ //spe3 
				if((status_of_submitted_card[0]==1)&&(status_of_submitted_card[1]==1)&&(status_of_submitted_card[2]==0x0001)){
					error=0;
					now_jocker_single=0;
					flash=1;
					ps[now_player].spe3++;
					if(WINDOW_TYPE==6 && (flag_wait_type!=3)){
						stage_card[0][7]+=100;
						initialize_window2(stage_card, old_stage_card, players_card,&win_info,now_pass, player_name, sekijun, &accept_flag, mibun);
						stage_card[0][7]-=100;
						wait_control(&win_info,&g_win_info, &flag_wait_type, &graph_mode, &g_flag);
					}
				}else{
					error|=1;
				}
			}
			if(debug){printf("error04 = %i: jocker \n",error);} //DEBUG

			if(error==0){
				switch (status_of_submitted_card[3]){
					case 0: // pass
						break;
					case 1: // single
						break;
					case 2: // pair
						break;
					case 3: // kaidan
						if(RULE_KAIDAN){
							break;
						}else{
							error|=1;
						}
						break;
					case 4: // jocker only
						now_jocker_single=1;
						now_shibari=0;
						break;
				}
			}
			if(debug){printf("error05 = %i: kaidan \n",error);} //DEBUG
		
		
			if(error || (status_of_submitted_card[1]==0)){
				if(error){
					if(debug){printf("error 01 - error = %i - number of card is %i\n",error,status_of_submitted_card[1]);}
					i=8;
				}else{
					if(debug){printf("OK 02\n");}
					i=9;
				}
				j=0; while(j<=0){
					j=tn_int_write(client_sockfd[now_player],i,protocol_version);
					if(debug){printf("ack roop %i\n",j);}
				};
				if(debug){printf("ack to player %i -- number %i\n",now_player,i);} //DEBUG
				now_pass[now_player]=1;
				number_renzoku_pass++;
			}else{
				if(debug){printf("OK 03\n");}
				bcopy(status_of_submitted_card,status_of_stages_card,sizeof(status_of_submitted_card));
				i=9;
				j=0;while(j<=0){j=tn_int_write(client_sockfd[now_player],i,protocol_version);if(debug){printf("ack roop %i",j);}};
				if(debug){printf("ack to player %i -- number %i\n",now_player,i);} //DEBUG
				bcopy(stage_card,old_stage_card,sizeof(work_card));
				bcopy(work_card,stage_card,sizeof(work_card));
				drop_card_flag(players_card[now_player],work_card);
				last_player=now_player;
				accept_flag=last_player;
				number_renzoku_pass=0;

				count_card_num(players_card[now_player], &number_of_card);
				if((number_of_card==0)&&(goal_flag[now_player]==0)){
					j=0;
					for(i=0;i<=4;i++){
						count_card_num(players_card[i], &number_of_card);
						j=j+(number_of_card==0);
					}
					goal_flag[now_player]=1;
					mibun[j-1]=now_player;
					now_number_of_goal=j;
				}

				if(now_jocker_single==0){
					if(RULE_KAKUMEI){
						if((status_of_submitted_card[3]==2)&&(status_of_submitted_card[1]>=4)){
							now_kakumei=(now_kakumei+1)%2;
							ps[now_player].kakumei++;
							
							if((WINDOW_TYPE==6) &&(flag_wait_type!=3)){
								stage_card[5][6]+=100;
								initialize_window2(stage_card, old_stage_card, players_card,&win_info,now_pass, player_name, sekijun, &accept_flag, mibun);
								stage_card[5][6]-=100;
								wait_control(&win_info,&g_win_info, &flag_wait_type, &graph_mode, &g_flag);
							}
						}
						if((status_of_submitted_card[3]==3)&&(status_of_submitted_card[1]>=5)){
							now_kakumei=(now_kakumei+1)%2;
							ps[now_player].kakumei++;
						}
					}
					if(RULE_8GIRI){
						if(check_special_card(8-2,status_of_submitted_card,0)){
							flash=1;
							ps[now_player].eightGiriCnt++;
							ps[now_player].eightGiriTurnSum+=count_turn;
						}
					}
					if(RULE_6REVERS){
						if(check_special_card(6-2,status_of_submitted_card,0)){
							if(now_muki==1){
								now_muki=4;
							}else{
								now_muki=1;
							}
						}
					}
					if(RULE_5TOBI){
						if(check_special_card(5-2,status_of_submitted_card,0)){
							i=get_seat(sekijun,now_player);
							now_player=sekijun[(i+now_muki)%5];
						}
					}
					if(RULE_11BACK){
						if(check_special_card(11-2,status_of_submitted_card,0)){
							now_11back=(now_11back+1)%2;	
						}
					}
				}

			}
	
			//convert status to table
			work_card[5][0] = 0;  
			work_card[5][1] = 0;
			work_card[5][2] = 0;
			work_card[5][3] = now_player; 
			work_card[5][4] = (status_of_stages_card[0]==-1);
			work_card[5][5] = now_11back; 
			work_card[5][6] = now_kakumei;
			work_card[5][7] = now_shibari;
			work_card[6][0] = count_card_num_r(players_card[0], &error);
			work_card[6][1] = count_card_num_r(players_card[1], &error);
			work_card[6][2] = count_card_num_r(players_card[2], &error);
			work_card[6][3] = count_card_num_r(players_card[3], &error);
			work_card[6][4] = count_card_num_r(players_card[4], &error);
			if(((now_number_of_games-1 )% FLASH_MIBUN_NUMBER)==0){
				for(i=0;i<=4;i++){
					work_card[6][i+5]=0;
				}
			}else{
				for(i=0;i<=4;i++){
					work_card[6][mibun[i]+5]=i;
				}
			}
			for(i=10;i<=14;i++){
				work_card[6][i]=sekijun[i-10];
			}
			for(i=0;i<=4;i++){
				bcopy(work_card[5], stage_card[5], 2*sizeof(work_card[5]));
			}
			bcopy(work_card[5], stage_card[5], 2*sizeof(work_card[5]));
			players_card[now_player][5][2]=1;

			// send information "stage cards" to clients.
			tn_card_write(client_sockfd[0],stage_card,protocol_version);
			tn_card_write(client_sockfd[1],stage_card,protocol_version);
			tn_card_write(client_sockfd[2],stage_card,protocol_version);
			tn_card_write(client_sockfd[3],stage_card,protocol_version);
			tn_card_write(client_sockfd[4],stage_card,protocol_version);

			if(number_renzoku_pass>=20){
				if(debug){printf("renzoku pass \n");}
				// srand((unsigned)time(NULL));	
				while(now_number_of_goal<=3){
					j=(int)(tn_rand_gen(RAND_TYPE)*(5-now_number_of_goal)+1);
					i=0;
					k=0;
					while(i<j){
						k++;
						if(goal_flag[k-1]==0){
							i++;
						}
					}
					goal_flag[k-1]=1;
					mibun[now_number_of_goal]=k-1;
					now_number_of_goal++;
				}
			}


			if((now_pass[0]+now_pass[1]+now_pass[2]+now_pass[3]+now_pass[4])>=(5-now_number_of_goal)){
				flash=1;	
			}

			if((WINDOW_TYPE==6)&&(flag_wait_type!=3)){
				initialize_window2(stage_card, old_stage_card, players_card,&win_info,now_pass, player_name, sekijun, &accept_flag, mibun);
				wait_control(&win_info,&g_win_info, &flag_wait_type, &graph_mode, &g_flag);
			}

			if(flash){
				flash=0;	
				now_11back=0;	
				now_player=last_player;
				now_jocker_single=0;
				now_shibari=0;   // flag of jocker
				bzero(now_pass,sizeof(now_pass));
				status_of_stages_card[0]=-1;
				status_of_stages_card[1]=-1;
				status_of_stages_card[2]=-1;
				status_of_stages_card[3]=-1;
				bcopy(stage_card,old_stage_card,sizeof(stage_card));
				bzero(stage_card,sizeof(stage_card));

				count_card_num(players_card[now_player], &number_of_card);
				while(number_of_card==0){
					if(debug){printf("now_player search %i \n",now_player);} //DEBUG
					i=get_seat(sekijun,now_player);
					now_player=sekijun[(i+now_muki)%5];
					count_card_num(players_card[now_player], &number_of_card);
				}
				last_player=now_player;
				if(debug){printf("flashed 01==>next player is %i %\n",now_player);}
			}else{
				if(debug){printf("no flash==>%i %i %i \n",now_player, now_muki,(now_player+now_muki)%5);}
				i=get_seat(sekijun,now_player);
				now_player=sekijun[(i+now_muki)%5];

				count_card_num(players_card[now_player], &number_of_card);
				while(((number_of_card==0)||(now_pass[now_player]==1))){
					if(debug){printf("now_player search %i \n",now_player);} //DEBUG
					i=get_seat(sekijun,now_player);
					now_player=sekijun[(i+now_muki)%5];
					count_card_num(players_card[now_player], &number_of_card);
				}
			}
			if(debug){printf("game is contineous = %i \n",now_number_of_goal);}

			if(now_number_of_goal==4){
				accept_flag=10;
				if(now_number_of_games==GAME_NUMBER){ // send a information "all game is overd" to clients.
					i=2;
					for(j=0;j<=4;j++){
						tn_int_write(client_sockfd[j],i,protocol_version);
					}
				}else{  // send a information "One game is overd" to clients.
					i=1;
					for(j=0;j<=4;j++){
						tn_int_write(client_sockfd[j],i,protocol_version);
					}
				}
			}else{
				i=0;
				for(j=0;j<=4;j++){
					tn_int_write(client_sockfd[j],i,protocol_version);
				}
			}
			sum_of_turn++;
		}// elihw of 1 game
		for(i=0;i<=4;i++){
			if(goal_flag[i]==0){
				mibun[4]=i;
			}
		}
		if(debug){printf("end of 1 game\n");}
		if(WINDOW_TYPE==2){
			printf("================ game %i \n",now_number_of_games);
			printf("daihugou   %i\n",mibun[0]);
			printf("hugou      %i\n",mibun[1]);
			printf("heimin     %i\n",mibun[2]);
			printf("himmin     %i\n",mibun[3]);
			printf("daihinnmin %i\n",mibun[4]);
			printf("--------- total point\n");
		}
		for(i=0;i<=4;i++){
			human[mibun[i]]=i;
			point[mibun[i]][now_number_of_games]=point[mibun[i]][now_number_of_games-1]+(5-i);
		}
		if(WINDOW_TYPE==2){
			printf("%u %u %u %u %u\n",point[0][now_number_of_games],point[1][now_number_of_games],point[2][now_number_of_games],point[3][now_number_of_games],point[4][now_number_of_games]);
		}
		fprintf(fp2,"%i %i %i %i %i\n",human[0],human[1],human[2],human[3],human[4]);
		fprintf(fp,"%u %u %u %u %u\n",point[0][now_number_of_games],point[1][now_number_of_games],point[2][now_number_of_games],point[3][now_number_of_games],point[4][now_number_of_games]);

		if(GRAPH_WINDOW==1 || GRAPH_WINDOW==2){
			switch(graph_mode){
				case 0: graph_test(point,now_number_of_games,&g_win_info, player_name, &g_flag);break;
				case 1: graph_test2(point,now_number_of_games,&win_info, &g_win_info, player_name, &g_flag);break;
			}
			XFlush( g_win_info.display );
		}//fi GRAPH_WINDOW
		game_count++;

	} // rof now_number_of_games

	// player statistics
	printf("Player Statistics\n");
	for (i=0;i<=4;i++) {
		printf("%s\n", player_name[i]);
		if (ps[i].getStage !=0) {	
		printf("average of card strength:\t\t%f\n", (double)ps[i].cardStrength/(double)ps[i].getStage);
		printf("average of shibari:\t\t%f\n", (double)ps[i].shibari/(double)ps[i].shibariCnt);
		printf("average of fukusuu(done/chance):\t\t%f\n", (double)ps[i].fukusuu/(double)ps[i].getStage);
		printf("average of kaidan(done/chance):\t\t%f\n", (double)ps[i].kaidan/(double)ps[i].getStage);
		}	
		if (ps[i].jokerCnt != 0) 
			printf("average of joker turn:\t\t%d\n", ps[i].jokerTurnSum/ps[i].jokerCnt);
		if (ps[i].eightGiriCnt != 0)
			printf("average of 8giri turn:\t\t%d\n", ps[i].eightGiriTurnSum/ps[i].eightGiriCnt);
		printf("average of spe3:\t\t%d\n", ps[i].spe3);
		printf("average of kakumei:\t\t%d\n",ps[i].kakumei);
	}
	printf("\n");
	printf("all turn ;\t\t%d , all game: \t\t%d\n",sum_of_turn, game_count);
	printf("average of turn:\t\t%f\n",(double) sum_of_turn / (double) game_count);

	/*************/ 
	/* game over */ 
	/*************/ 
	printf("All games are overed \n");
	for(i=0;i<=4;i++){	
		shutdown(client_sockfd[i], 2);
		close(client_sockfd[i]);
	}

	if(WINDOW_TYPE==0 || WINDOW_TYPE==1 || WINDOW_TYPE==3 || WINDOW_TYPE==4 || WINDOW_TYPE==5 || WINDOW_TYPE==6){ 
		stop_control(&win_info);
	}

	shutdown(sockfd, 2);
	close(sockfd);
	fclose(fp);
	fclose(fp2);

}//niam
