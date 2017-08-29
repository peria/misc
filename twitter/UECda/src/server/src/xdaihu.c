/*
  Grover's amplitude
 */
#include<sys/time.h>
#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include<time.h>
#include<string.h>
#include<ctype.h>
#include"xdaihu.h"
#include"external.h"

int isint(char* a){
	/************************************************************/
	/* input  : A string                                        */
	/* return : If a input string are integer and               */
	/*        : the length of the string smaller than 5 then 1, */
	/*        : otherwize 0.                                    */
	/* destroy: none                                            */
	/************************************************************/
	int i;

	if(strlen(a)>=6){
		return 0;
	}
	for(i=0;i<strlen(a);i++){
		if(isdigit(a[i])==0){
			return 0;
		}
	}
	return 1;
}

char* num_to_str(char* str,int flag){
	/************************************************************/
	/* input  : str : strings                                   */
	/*        : flag : an integer which corresponding to str.   */
	/* return : an strings corresponding to str and flag        */
	/* destroy: none                                            */
	/************************************************************/
	if(strcmp(str,"WINDOW_TYPE")==0){	
		switch(flag){ 
			case 0:
				return "SMALL";
			case 1:
				return "BIG";
			case 2:
				return "CONSOLE";
			case 3:
				return "MIDDLE_EXTEND";
			case 4:
				return "MIDDLE";
			case 5:
				return "MIDDLE_EXTEND";
			case 6:
				return "MIDDLE_XGA";
		}
	}
	return "ERROR";
}

char* one_to_yes(int flag){
	/************************************************************/
	/* input  : flag : integer                                  */
	/* return : If flag==1 the retunr YES otherwize NO          */
	/* destroy: none                                            */
	/************************************************************/
	if(flag){ 
		return "YES";
	}else{
		return "NO";
	}
}

void strupr2(char *dat){
	while(*dat != '\0'){ 
		//if(*dat >= 'a' && *dat <= 'z')
		// 	*dat = *dat - 0x20;
		*dat=toupper(*dat);
		dat++; 
	}
}


unsigned long MyColor(Display *display, char *color){
	Colormap cmap;
	XColor c0, c1;

	cmap = DefaultColormap(display,0);

	XAllocNamedColor(display, cmap,color,&c1, &c0);
	return (c1.pixel);
}


void graph_test(unsigned int point[5][15000],int game_num,g_window_info *g_win_info, char player_name[5][15], int *flag){

	int i,j,k;
	unsigned long  color[5];
	char value[6][10];
	int font_offset=10;
	int y_font_offset=20;
	double average[5];

	color[0]=MyColor(g_win_info->display,"rgb:99/00/ff");
	color[1]=MyColor(g_win_info->display,"red");
	color[2]=MyColor(g_win_info->display,"rgb:00/00/ff");
	color[3]=MyColor(g_win_info->display,"brown");
	color[4]=MyColor(g_win_info->display,"rgb:00/75/00");

	for(i=0;i<=4;i++){
		if(g_win_info->y_border<=g_win_info->y_scale*point[i][game_num]){
			g_win_info->y_scale=g_win_info->y_scale/2;
			*flag=1;
		}
		if(g_win_info->x_border<=g_win_info->x_scale*game_num){
			g_win_info->x_scale=g_win_info->x_scale/2;
			*flag=1;
		}
	}
	if(game_num==1){
		*flag=1;
		//printf("debug: graph 1st\n");
	}

	if(*flag){
		sprintf(value[0],"%.0fpt \0",(g_win_info->y_border)/g_win_info->y_scale/2);
		sprintf(value[1],"%.0fpt \0",(g_win_info->y_border)/g_win_info->y_scale);
		sprintf(value[2],"%.0f\0",1000/g_win_info->x_scale);
		sprintf(value[3],"%.0f\0",750/g_win_info->x_scale);
		sprintf(value[4],"%.0f\0",500/g_win_info->x_scale);
		sprintf(value[5],"%.0f\0",250/g_win_info->x_scale);

		XSetForeground(g_win_info->display,g_win_info->gc,g_win_info->background); // clean pixmap
		XSetBackground(g_win_info->display,g_win_info->gc,g_win_info->background);
		XFillRectangle(g_win_info->display, g_win_info->pixmap, g_win_info->gc, 0,0, g_win_info->window_width, g_win_info->window_height);

		XSetForeground(g_win_info->display,g_win_info->gc,g_win_info->foreground);
		XSetBackground(g_win_info->display,g_win_info->gc,g_win_info->background);

		for(i=1;i<=4;i++){ // draw axis 
			XDrawLine(g_win_info->display, g_win_info->pixmap, g_win_info->gc, i*250,0,
				i*250, g_win_info->window_height-20);
			XDrawImageString(g_win_info->display, g_win_info->pixmap, g_win_info->gc, 
				1000-250*(i-1)-font_offset, g_win_info->window_height, value[i+1], strlen(value[i+1]));
		}
		for(i=1;i<=2;i++){
			XDrawLine(g_win_info->display, g_win_info->pixmap, g_win_info->gc, 0,g_win_info->window_height-g_win_info->y_border*i/2,
				g_win_info->window_width, g_win_info->window_height-g_win_info->y_border*i/2);
			XDrawImageString(g_win_info->display, g_win_info->pixmap, g_win_info->gc, 
				0, g_win_info->window_height-g_win_info->y_border*i/2+font_offset, value[i-1], strlen(value[i-1]));

		}

		XSetForeground(g_win_info->display,g_win_info->gc,g_win_info->foreground); // write player number
		for(i=0;i<=4;i++){
			sprintf(value[0],":%6u \0",point[i][game_num]);
			XDrawImageString(g_win_info->display, g_win_info->pixmap, g_win_info->gc, 
				g_win_info->window_width-280, g_win_info->window_height-180+30*i+font_offset, player_name[i], strlen(player_name[i]));
			XDrawImageString(g_win_info->display, g_win_info->pixmap, g_win_info->gc, 
				g_win_info->window_width-180, g_win_info->window_height-180+30*i+font_offset, value[0], strlen(value[0]));
		}

		for(k=0;k<=4;k++){ // old graph
			XSetForeground(g_win_info->display,g_win_info->gc,color[k]);
			for(j=1;j<=game_num;j++){ 
				XDrawLine(g_win_info->display, g_win_info->pixmap, g_win_info->gc, 
					g_win_info->x_scale*(j-1),g_win_info->window_height-g_win_info->y_scale*point[k][j-1],
					g_win_info->x_scale*j, g_win_info->window_height-g_win_info->y_scale*point[k][j]);
			}
			// player name's color
			XFillRectangle(g_win_info->display, g_win_info->pixmap, g_win_info->gc, 
				g_win_info->window_width-300,g_win_info->window_height-195+30*k, 10, 30);
			XFillRectangle(g_win_info->display, g_win_info->pixmap, g_win_info->gc, 
				g_win_info->window_width-300,g_win_info->window_height-195+30*(k+1)-2,
				120,2);
		}


	}


	for(i=0;i<=4;i++){ // draw line
		XSetForeground(g_win_info->display,g_win_info->gc,color[i]);
		XDrawLine(g_win_info->display, g_win_info->pixmap, g_win_info->gc, g_win_info->x_scale*(game_num-1),g_win_info->window_height-g_win_info->y_scale*point[i][game_num-1],
			g_win_info->x_scale*game_num, g_win_info->window_height-g_win_info->y_scale*point[i][game_num]);

		XSetForeground(g_win_info->display,g_win_info->gc,g_win_info->foreground); // write point
		sprintf(value[0],":%6u  \0",point[i][game_num]);
		XDrawImageString(g_win_info->display, g_win_info->pixmap, g_win_info->gc, 
			//920, g_win_info->window_height-180+30*i+font_offset, value[0], strlen(value[0]));
			g_win_info->window_width-180, g_win_info->window_height-180+30*i+font_offset, value[0], strlen(value[0]));
	}

	if((!(game_num%5) || *flag)&&(game_num!=1)){
		for(i=0;i<=4;i++){
			average[i]=((double)point[i][game_num]-(double)point[i][game_num-10])/(double)10;
			sprintf(value[0],"%5.1f\0",average[i]);
			XDrawImageString(g_win_info->display, g_win_info->pixmap, g_win_info->gc, 
				//1030, g_win_info->window_height-180+30*i+font_offset, value[0], strlen(value[0]));
				g_win_info->window_width-70, g_win_info->window_height-180+30*i+font_offset, value[0], strlen(value[0]));
		}
	}
	*flag=0;
	XCopyArea(g_win_info->display, g_win_info->pixmap, g_win_info->window, g_win_info->gc,  // Copy pixmap to window
		0,0, g_win_info->window_width, g_win_info->window_height, 0,0); 
}

void graph_initialize(g_window_info *g_win_info, char player_name[5][15]){

	int i;
	char value[6][10];
	unsigned long  color[5];
	int font_offset=10;
	int y_font_offset=20;
	
	color[0]=MyColor(g_win_info->display,"rgb:99/00/ff");
	color[1]=MyColor(g_win_info->display,"red");
	color[2]=MyColor(g_win_info->display,"rgb:00/00/ff");
	color[3]=MyColor(g_win_info->display,"brown");
	color[4]=MyColor(g_win_info->display,"rgb:00/75/00");

	sprintf(value[0],"%.0fpt \0",(g_win_info->y_scale*50)/g_win_info->y_scale);
	sprintf(value[1],"%.0fpt \0",(g_win_info->y_scale*100)/g_win_info->y_scale);
	sprintf(value[2],"%.0f\0",1000/g_win_info->x_scale);
	sprintf(value[3],"%.0f\0",750/g_win_info->x_scale);
	sprintf(value[4],"%.0f\0",500/g_win_info->x_scale);
	sprintf(value[5],"%.0f\0",250/g_win_info->x_scale);

	XSetForeground(g_win_info->display,g_win_info->gc,g_win_info->background); // clean pixmap
	XSetBackground(g_win_info->display,g_win_info->gc,g_win_info->background);
	XFillRectangle(g_win_info->display, g_win_info->pixmap, g_win_info->gc, 0,0, g_win_info->window_width, g_win_info->window_height);

	XSetForeground(g_win_info->display,g_win_info->gc,g_win_info->foreground);
	XSetBackground(g_win_info->display,g_win_info->gc,g_win_info->background);

	for(i=1;i<=4;i++){ // draw axis 
		XDrawLine(g_win_info->display, g_win_info->pixmap, g_win_info->gc, i*250,0,
			i*250, g_win_info->window_height-20);
		XDrawImageString(g_win_info->display, g_win_info->pixmap, g_win_info->gc, 
			1000-250*(i-1)-font_offset, g_win_info->window_height, value[i+1], strlen(value[i+1]));
	}
	for(i=1;i<=2;i++){
		XDrawLine(g_win_info->display, g_win_info->pixmap, g_win_info->gc, 0,g_win_info->window_height-(g_win_info->y_scale*50)*i,
			g_win_info->window_width, g_win_info->window_height-(g_win_info->y_scale*50)*i);
		XDrawImageString(g_win_info->display, g_win_info->pixmap, g_win_info->gc, 
			0, g_win_info->window_height-(g_win_info->y_scale*50)*i+font_offset, value[i-1], strlen(value[i-1]));
	}

	XSetForeground(g_win_info->display,g_win_info->gc,g_win_info->foreground); // write player number
	for(i=0;i<=4;i++){
		XDrawImageString(g_win_info->display, g_win_info->pixmap, g_win_info->gc, 
			g_win_info->window_width-280, g_win_info->window_height-180+30*i+font_offset, player_name[i], strlen(player_name[i]));
	}

	for(i=0;i<=4;i++){
		XSetForeground(g_win_info->display,g_win_info->gc,color[i]);
		XFillRectangle(g_win_info->display, g_win_info->pixmap, g_win_info->gc, 
			g_win_info->window_width-300, g_win_info->window_height-195+30*i, 10, 30);
		XFillRectangle(g_win_info->display, g_win_info->pixmap, g_win_info->gc, 
			g_win_info->window_width-300,g_win_info->window_height-195+30*(i+1)-2,
			120,2);
	}
	
	XCopyArea(g_win_info->display, g_win_info->pixmap, g_win_info->window, g_win_info->gc,  // Copy pixmap to window
		0,0, g_win_info->window_width, g_win_info->window_height, 0,0); 
}

void initialize_window(int stage_card[8][15], int old_stage_card[8][15], int players_card[5][8][15], c_window_info *win_info,int now_pass[5], char player_name[5][15], int sekijun[5]){

	int i,j,k;
	int xflag,yflag;
	int number_of_stage_card;
	int number_of_player_card;
	Colormap cmap;
	XColor color, exact;
	unsigned long green, blue, red;


	// initialize stage on piro

	if(stage_card[5][6]){  // kakumei 
		XCopyArea(win_info->display, win_info->screen4, win_info->screen0, win_info->gc, 0,0,win_info->window_width, win_info->window_height, 0,0); // delete old card on stage
	}else{
		XCopyArea(win_info->display, win_info->screen1, win_info->screen0, win_info->gc, 0,0,win_info->window_width, win_info->window_height, 0,0); // delete old card on stage
	}
	if(stage_card[5][7]){  // shibari
		XCopyArea(win_info->display,win_info->mibun_win,win_info->screen0,win_info->gc,0,win_info->card_height*7,win_info->card_width,win_info->card_height,win_info->on_stage_card_x_offset+win_info->card_width*10,win_info->on_stage_card_y_offset-(win_info->window_type==3)*10);
	}
	for(i=0;i<=4;i++){ // mibun and pass and turn
		XCopyArea(win_info->display,win_info->mibun_win,win_info->screen0,win_info->gc,0,win_info->card_height*6,win_info->card_width,win_info->card_height,5+win_info->card_width+win_info->on_stage_player_x_offset,win_info->on_stage_player_size*get_seat(sekijun,players_card[0][5][3])+win_info->on_stage_player_y_offset); // turn	
		if(now_pass[i]==1){
			XCopyArea(win_info->display,win_info->mibun_win,win_info->screen0,win_info->gc,0,0,win_info->card_width,win_info->card_height,win_info->on_stage_card_x_offset+win_info->card_width+win_info->on_stage_player_x_offset,win_info->on_stage_player_size*get_seat(sekijun,i)+win_info->on_stage_player_y_offset); //pass	
		}
	}

	if(win_info->window_type==3 || win_info->window_type==5){
		// initialize stage for old cards 
		number_of_stage_card=0;
		for(j=0;j<=14;j++){
			for(i=0;i<=4;i++){
				switch (old_stage_card[i][j]){
					case 1:
						xflag=(j-1)*win_info->card_width;
						yflag=i*win_info->card_height;
						XCopyArea(win_info->display, win_info->screen2, win_info->screen0, win_info->gc, xflag,yflag,win_info->card_width,win_info->card_height,win_info->on_stage_card_x_offset+win_info->card_width*number_of_stage_card+10,win_info->on_stage_card_y_offset-20);
						number_of_stage_card++;
						break;
					case 2:
						xflag=0*win_info->card_width;
						yflag=4*win_info->card_height;
						XCopyArea(win_info->display, win_info->screen2, win_info->screen0, win_info->gc, xflag,yflag,win_info->card_width,win_info->card_height,win_info->on_stage_card_x_offset+win_info->card_width*number_of_stage_card+10,win_info->on_stage_card_y_offset-20);
						number_of_stage_card++;
				}
			}
		}
	}

	if(win_info->window_type==5){
		for(i=0;i<=4;i++){
			XDrawImageString(win_info->display, win_info->screen0, win_info->gc,
				5, 220+115*get_seat(sekijun,i), player_name[i], strlen(player_name[i]));
		}
	}


	// initialize stage on piro
	number_of_stage_card=0;
	for(j=0;j<=14;j++){
		for(i=0;i<=4;i++){
			switch (stage_card[i][j]){
				case 1:
					xflag=(j-1)*win_info->card_width;
					yflag=i*win_info->card_height;
					XCopyArea(win_info->display, win_info->screen2, win_info->screen0, win_info->gc, xflag,yflag,win_info->card_width,win_info->card_height,win_info->on_stage_card_x_offset+win_info->card_width*number_of_stage_card,win_info->on_stage_card_y_offset);
					number_of_stage_card++;
					break;
				case 2:
					xflag=0*win_info->card_width;
					yflag=4*win_info->card_height;
					XCopyArea(win_info->display, win_info->screen2, win_info->screen0, win_info->gc, xflag,yflag,win_info->card_width,win_info->card_height,win_info->on_stage_card_x_offset+win_info->card_width*number_of_stage_card,win_info->on_stage_card_y_offset);
					number_of_stage_card++;
			}
		}
	}
	if(win_info->window_type==3 || win_info->window_type==5){
		if(number_of_stage_card==0){
			XCopyArea(win_info->display, win_info->screen2, win_info->screen0, win_info->gc, 
				5*win_info->card_width, 4*win_info->card_height,
				4*win_info->card_width,win_info->card_height,
				win_info->on_stage_card_x_offset+win_info->card_width*4,win_info->on_stage_card_y_offset-20);
		}
	}

	// initializ cards on players cards
	for(k=0;k<=4;k++){
		number_of_player_card=0;
		for(j=1;j<=13;j++){
			for(i=0;i<=3;i++){
				if(players_card[k][i][j]==1){
					xflag=(j-1)*win_info->card_width;
					yflag=i*win_info->card_height;
					XCopyArea(win_info->display, win_info->screen2, win_info->screen0, win_info->gc, xflag,yflag,win_info->card_width,win_info->card_height,win_info->card_offset_x+win_info->card_width*number_of_player_card+win_info->on_stage_player_x_offset,win_info->on_stage_player_size*get_seat(sekijun,k)+win_info->on_stage_player_y_offset);
					number_of_player_card++;
				}
			}
		}
		if(players_card[k][4][1]==2){
			xflag=0*win_info->card_width;
			yflag=4*win_info->card_height;
			XCopyArea(win_info->display, win_info->screen2, win_info->screen0, win_info->gc, xflag,yflag,win_info->card_width,win_info->card_height,win_info->card_offset_x+win_info->card_width*number_of_player_card+win_info->on_stage_player_x_offset,win_info->on_stage_player_size*get_seat(sekijun,k)+win_info->on_stage_player_y_offset);
			number_of_player_card++;
		}
	}
	XCopyArea(win_info->display, win_info->screen0, win_info->screen3, win_info->gc, 0,0,win_info->window_width, win_info->window_height, 0,0); // delete old card on stage
	XFlush( win_info->display );

}

void wait_control(c_window_info *win_info, g_window_info *g_win_info, int *flag_wait_type, int *graph_mode, int *g_flag){

	XEvent event;
	int i;
	int press;
	struct timeval waitval;

	while( (XEventsQueued( win_info->display, QueuedAfterFlush ) !=0) ){ // flash XEvent que
		XNextEvent(win_info->display, &event);
//		XCheckMaskEvent(win_info->display,  ButtonPressMask,&event);
		switch(event.type){
			case ButtonPress :
				press=event.xbutton.y/50;
				if((0<=press) && (press <=3)){ 
					*flag_wait_type=press;
				}else if(press==4){
					*graph_mode=(*graph_mode+1)%2;
					*g_flag=1;
				}
				break;
		}
	}

	switch (*flag_wait_type){
		case 0:
			XRaiseWindow(win_info->display,win_info->screen3);
			XLowerWindow(win_info->display,g_win_info->window);
			break;
		case 1:
			waitval.tv_sec  = 1;
			waitval.tv_usec = 00;
			select(0,NULL,NULL,NULL,&waitval);
			XRaiseWindow(win_info->display,win_info->screen3);
			XLowerWindow(win_info->display,g_win_info->window);
			break;
		case 2: // wait press bottun
			while( XEventsQueued( win_info->display, QueuedAfterFlush ) != 0){ // flash XEvent que
				XNextEvent(win_info->display, &event);
			}
			XRaiseWindow(win_info->display,win_info->screen3);
			XLowerWindow(win_info->display,g_win_info->window);
			XMaskEvent(win_info->display,  ButtonPressMask,&event);
			*flag_wait_type=event.xbutton.y/50;
			break;
		case 3: // XRaiseWindow(win_info->display,g_win_info->window);
			XLowerWindow(win_info->display,win_info->screen3);
			break;
	}
}

void stop_control(c_window_info *win_info){
	XEvent event;
	int i;

	XMaskEvent(win_info->display,  ButtonPressMask,&event);
	XMaskEvent(win_info->display,  ButtonPressMask,&event);
	XMaskEvent(win_info->display,  ButtonPressMask,&event);
}

void graph_test2(unsigned int point[5][15000],int game_num,c_window_info *win_info, g_window_info *g_win_info, char player_name[5][15], int *flag){
	int i,j,k;
	unsigned long  color[5];
	char value[6][10];
	int font_offset=10;
	int y_font_offset=20;
	double average[5];

	color[0]=MyColor(g_win_info->display,"rgb:99/00/ff");
	color[1]=MyColor(g_win_info->display,"red");
	color[2]=MyColor(g_win_info->display,"rgb:00/00/ff");
	color[3]=MyColor(g_win_info->display,"brown");
	color[4]=MyColor(g_win_info->display,"rgb:00/75/00");

	for(i=0;i<=4;i++){
		if(g_win_info->y_border<=g_win_info->y_scale2*point[i][game_num]){
			g_win_info->y_scale2=g_win_info->y_scale2/2;
			*flag=1;
		}
	}
	if(game_num==1){
		*flag=1;
	}

	if(*flag){
		sprintf(value[0],"%.0f\0",0*g_win_info->y_scale2);
		sprintf(value[1],"%.0f\0",(g_win_info->y_border)/g_win_info->y_scale2/2);
		sprintf(value[2],"%.0f\0",(g_win_info->y_border)/g_win_info->y_scale2);

		XSetForeground(g_win_info->display,g_win_info->gc,g_win_info->background); // clean pixmap
		XSetBackground(g_win_info->display,g_win_info->gc,g_win_info->background);
		XFillRectangle(g_win_info->display, g_win_info->pixmap2, g_win_info->gc, 0,0, g_win_info->window_width, g_win_info->window_height);

		XSetForeground(g_win_info->display,g_win_info->gc,g_win_info->foreground);
		XSetBackground(g_win_info->display,g_win_info->gc,g_win_info->background);

		for(i=0;i<=2;i++){ // draw axis 
			XDrawLine(g_win_info->display, g_win_info->pixmap2, g_win_info->gc, 150+(g_win_info->y_border/2*i),0,
				150+(g_win_info->y_border/2*i), g_win_info->window_height-20);
			XDrawImageString(g_win_info->display, g_win_info->pixmap2, g_win_info->gc, 
				150+(g_win_info->y_border/2*i)-font_offset, g_win_info->window_height, value[i], strlen(value[i]));
		}

		XSetForeground(g_win_info->display,g_win_info->gc,g_win_info->foreground); // write player number
		for(i=0;i<=4;i++){
			sprintf(value[0],":%6u \0",point[i][game_num]);
			XDrawImageString(g_win_info->display, g_win_info->pixmap2, g_win_info->gc, 
				30, 125+130*i, player_name[i], strlen(player_name[i]));
			XCopyArea(g_win_info->display, win_info->screen6, g_win_info->pixmap2, g_win_info->gc,  // Copy pixmap to window
				650+125*i,710, 125, 115, 15,10+130*i); 
		}
	}

	for(i=0;i<=4;i++){ // draw line
		XSetForeground(g_win_info->display,g_win_info->gc,color[i]);
		//XFillRectangle(g_win_info->display, g_win_info->pixmap2, g_win_info->gc, 15,650+130*i, 130, 105 );
		XFillRectangle(g_win_info->display, g_win_info->pixmap2, g_win_info->gc, 150,20+130*i, g_win_info->y_scale2*point[i][game_num], 100 );
		XSetForeground(g_win_info->display,g_win_info->gc,g_win_info->foreground); // write point
		sprintf(value[0],"%6u  \0",point[i][game_num]);
		XDrawImageString(g_win_info->display, g_win_info->pixmap2, g_win_info->gc, 
			g_win_info->window_width-100, 80+130*i, value[0], strlen(value[0]));
	}
	if((!(game_num%5) || *flag)&&(game_num!=1)){
		for(i=0;i<=4;i++){
			average[i]=((double)point[i][game_num]-(double)point[i][game_num-10])/(double)10;
			sprintf(value[0],"%5.1f\0",average[i]);
			XDrawImageString(g_win_info->display, g_win_info->pixmap2, g_win_info->gc, 
				g_win_info->window_width-100, 110+130*i, value[0], strlen(value[0]));
		}
	}
	*flag=0;
	XCopyArea(g_win_info->display, g_win_info->pixmap2, g_win_info->window, g_win_info->gc,  // Copy pixmap to window
		0,0, g_win_info->window_width, g_win_info->window_height, 0,0); 
}

void graph_initialize2(g_window_info *g_win_info, char player_name[5][15]){

	int i,j,k;
	int flag=1;
	unsigned long  color[5];
	char value[6][10];
	int font_offset=10;
	int y_font_offset=20;
	double average[5];

	color[0]=MyColor(g_win_info->display,"rgb:99/00/ff");
	color[1]=MyColor(g_win_info->display,"red");
	color[2]=MyColor(g_win_info->display,"rgb:00/00/ff");
	color[3]=MyColor(g_win_info->display,"brown");
	color[4]=MyColor(g_win_info->display,"rgb:00/75/00");

	for(i=0;i<=4;i++){
		if(g_win_info->y_border<=g_win_info->y_scale2*1){
			g_win_info->y_scale2=g_win_info->y_scale2/2;
			flag=1;
		}
	}

	if(flag){
		sprintf(value[0],"%.0f\0",0*g_win_info->y_scale2);
		sprintf(value[1],"%.0f\0",(g_win_info->y_border)/g_win_info->y_scale2/2);
		sprintf(value[2],"%.0f\0",(g_win_info->y_border)/g_win_info->y_scale2);

		XSetForeground(g_win_info->display,g_win_info->gc,g_win_info->background); // clean pixmap
		XSetBackground(g_win_info->display,g_win_info->gc,g_win_info->background);
		XFillRectangle(g_win_info->display, g_win_info->pixmap2, g_win_info->gc, 0,0, g_win_info->window_width, g_win_info->window_height);

		XSetForeground(g_win_info->display,g_win_info->gc,g_win_info->foreground);
		XSetBackground(g_win_info->display,g_win_info->gc,g_win_info->background);

		for(i=0;i<=2;i++){ // draw axis 
			XDrawLine(g_win_info->display, g_win_info->pixmap2, g_win_info->gc, 150+(g_win_info->y_border/2*i),0,
				150+(g_win_info->y_border/2*i), g_win_info->window_height-20);
			XDrawImageString(g_win_info->display, g_win_info->pixmap2, g_win_info->gc, 
				150+(g_win_info->y_border/2*i)-font_offset, g_win_info->window_height, value[i], strlen(value[i]));
		}

		XSetForeground(g_win_info->display,g_win_info->gc,g_win_info->foreground); // write player number

	}
	XCopyArea(g_win_info->display, g_win_info->pixmap2, g_win_info->window, g_win_info->gc,  // Copy pixmap to window
		0,0, g_win_info->window_width, g_win_info->window_height, 0,0); 

}

void initialize_window2(int stage_card[8][15], int old_stage_card[8][15], int players_card[5][8][15], c_window_info *win_info,int now_pass[5], char player_name[5][15], int sekijun[5], int *accept_flag, int mibun[5]){

	int i,j,k;
	int xflag,yflag;
	int number_of_stage_card;
	int tmp,tmp2,loser;
	int number_of_player_card_x;
	int number_of_player_card_y;
	int rightside;
	int sa;
	int seat_x[5],seat_y[5];
	unsigned long  color[5];
	int analyze_result[4];
	int analyze_result2[4];
	int joker_flag=0;
	Colormap cmap;
	unsigned long green, blue, red;

	// set table
	seat_x[0]=0;seat_y[0]=0;
	seat_x[1]=1;seat_y[1]=0;
	seat_x[2]=2;seat_y[2]=0;
	seat_x[3]=0;seat_y[3]=1;
	seat_x[4]=2;seat_y[4]=1;
	color[0]=MyColor(win_info->display,"rgb:99/00/ff");
	color[1]=MyColor(win_info->display,"red");
	color[2]=MyColor(win_info->display,"rgb:00/00/ff");
	color[3]=MyColor(win_info->display,"brown");
	color[4]=MyColor(win_info->display,"rgb:00/75/00");

	// analyze_stages_card
	analyze_card(stage_card,analyze_result,0);

	// initialize state lump

	if(stage_card[5][6]){  // kakumei 
		XCopyArea(win_info->display, win_info->screen4, win_info->screen0, win_info->gc, 0,0,900,700, 0,0); // delete old card on stage
		XCopyArea(win_info->display,win_info->screen5,win_info->screen0,win_info->gc,0,60,70,180,300+10+0,350+5);
	}else{
		XCopyArea(win_info->display, win_info->screen1, win_info->screen0, win_info->gc, 0,0,900,700, 0,0); // delete old card on stage
	}
	if(stage_card[5][0]){  //  card koukan
		XCopyArea(win_info->display,win_info->screen5,win_info->screen0,win_info->gc,680,60,70*4,180,300+10+0,350+5);
	}
	if(analyze_result[3]==2){ //pair
		XCopyArea(win_info->display,win_info->screen5,win_info->screen0,win_info->gc,140,60,70,180,300+10+70,350+5);
	}else if(analyze_result[3]==3){ // kaidan
		XCopyArea(win_info->display,win_info->screen5,win_info->screen0,win_info->gc,70,60,70,180,300+10+140,350+5);
	}
	if(stage_card[5][7]){  // shibari
		XCopyArea(win_info->display,win_info->screen5,win_info->screen0,win_info->gc,210,60,70,180,300+10+210,350+5);
	}
	XSetForeground(win_info->display,win_info->gc,color[players_card[0][5][3]]);
	XSetLineAttributes(win_info->display,win_info->gc,10,LineSolid,CapNotLast,JoinMiter);
	XDrawRectangle(win_info->display,win_info->screen0,win_info->gc,300*seat_x[get_seat(sekijun,players_card[0][5][3])], 350* seat_y[get_seat(sekijun,players_card[0][5][3])],300,350);
	XSetForeground(win_info->display,win_info->gc,win_info->foreground); // write point
	for(i=0;i<=4;i++){ // mibun and pass and turn
		if(now_pass[i]==1){
			XCopyArea(win_info->display,win_info->screen5,win_info->screen0,win_info->gc,750,0,150,90, 130+300*seat_x[get_seat(sekijun,i)], 5+350*seat_y[get_seat(sekijun,i)]); //pass	
		}
	}

	// name and figure
	XDrawImageString(win_info->display, win_info->screen0, win_info->gc, 160,90 , player_name[sekijun[0]], strlen(player_name[sekijun[0]]));
	XDrawImageString(win_info->display, win_info->screen0, win_info->gc, 460,90 , player_name[sekijun[1]], strlen(player_name[sekijun[1]]));
	XDrawImageString(win_info->display, win_info->screen0, win_info->gc, 760,90 , player_name[sekijun[2]], strlen(player_name[sekijun[2]]));
	XDrawImageString(win_info->display, win_info->screen0, win_info->gc, 160,440 , player_name[sekijun[3]], strlen(player_name[sekijun[3]]));
	XDrawImageString(win_info->display, win_info->screen0, win_info->gc, 760,440 , player_name[sekijun[4]], strlen(player_name[sekijun[4]]));
	XCopyArea(win_info->display,win_info->screen6,win_info->screen0,win_info->gc,130*sekijun[0],700,130,105, 10,5); // turn	
	XCopyArea(win_info->display,win_info->screen6,win_info->screen0,win_info->gc,130*sekijun[1],700,130,105, 310,5); // turn	
	XCopyArea(win_info->display,win_info->screen6,win_info->screen0,win_info->gc,130*sekijun[2],700,130,105, 610,5); // turn	
	XCopyArea(win_info->display,win_info->screen6,win_info->screen0,win_info->gc,130*sekijun[3],700,130,105, 10,355); // turn	
	XCopyArea(win_info->display,win_info->screen6,win_info->screen0,win_info->gc,130*sekijun[4],700,130,105, 610,355); // turn	

	// initialize stage on piro
		// old card
	analyze_card(old_stage_card,analyze_result2,0);
	if(analyze_result2[1]==0){
		XCopyArea(win_info->display, win_info->screen5, win_info->screen0, win_info->gc, 
			914,320,
			140,80,
			300+50,550);
	}else if(analyze_result2[1]<=3){
		rightside=300+50+(analyze_result2[1]-1)*70;
		sa=70;	
	}else{
		rightside=550-70;
		sa=(240-70)/(analyze_result2[1]-1);	
	};
	number_of_stage_card=0;
	for(j=14;j>=0;j--){
		for(i=4;i>=0;i--){
			switch (old_stage_card[i][j]){
				case 1:
					xflag=(j-1)*70+2;
					yflag=240+i*80 ;
					XCopyArea(win_info->display, win_info->screen5, win_info->screen0, win_info->gc, 
						xflag,yflag,
						70,80,
						rightside-number_of_stage_card*sa,550);
					number_of_stage_card++;
					break;
				case 2:
					xflag=910+2;
					yflag=240;
					XCopyArea(win_info->display, win_info->screen5, win_info->screen0, win_info->gc, 
						xflag,yflag,
						70,80,
						rightside-number_of_stage_card*sa,550);
					number_of_stage_card++;
			}
		}
	}

		// now card
	if(analyze_result[1]==0){
		XCopyArea(win_info->display, win_info->screen5, win_info->screen0, win_info->gc, 
			914,320,
			140,80,
			310,600);
	}else if(analyze_result[1]<=3){
		rightside=300+10+(analyze_result[1]-1)*70;
		sa=70;	
	}else{
		rightside=550-70;
		sa=(240-70)/(analyze_result[1]-1);	
	};
	number_of_stage_card=0;
	for(j=14;j>=0;j--){
		for(i=4;i>=0;i--){
			switch (stage_card[i][j]){
				case 1:
					xflag=(j-1)*70+2;
					yflag=240+i*80 ;
					XCopyArea(win_info->display, win_info->screen5, win_info->screen0, win_info->gc, 
						xflag,yflag,
						70,80,
						rightside-number_of_stage_card*sa,600);
					if((*accept_flag!=10)){
						XSetForeground(win_info->display,win_info->gc,color[players_card[0][5][3]]);
						XSetLineAttributes(win_info->display,win_info->gc,10,LineSolid,CapNotLast,JoinMiter);
						XDrawRectangle(win_info->display,win_info->screen0,win_info->gc,
							rightside-number_of_stage_card*sa+2,600+2,
							70-8,80-8);
						XSetForeground(win_info->display,win_info->gc,win_info->foreground); // write point
					}
					number_of_stage_card++;
					break;
				case 2:
					xflag=910+2;
					yflag=240;
					XCopyArea(win_info->display, win_info->screen5, win_info->screen0, win_info->gc, 
						xflag,yflag,
						70,80,
						rightside-number_of_stage_card*sa,600);
					if((*accept_flag!=10)){
						XSetForeground(win_info->display,win_info->gc,color[players_card[0][5][3]]);
						XSetLineAttributes(win_info->display,win_info->gc,10,LineSolid,CapNotLast,JoinMiter);
						XDrawRectangle(win_info->display,win_info->screen0,win_info->gc,
							rightside-number_of_stage_card*sa+2,600+2,
							70-8,80-8);
						XSetForeground(win_info->display,win_info->gc,win_info->foreground); // write point
					}
					joker_flag=1;
					number_of_stage_card++;
			}
		}
	}

	if(stage_card[0][7]>=100){ // spe3 tokubetu
		XCopyArea(win_info->display,win_info->screen5,win_info->screen0,win_info->gc,
			912,425,
			118,135,
			310,600+80-135);
		XSetForeground(win_info->display,win_info->gc,color[players_card[0][5][3]]);
		XSetLineAttributes(win_info->display,win_info->gc,10,LineSolid,CapNotLast,JoinMiter);
		XDrawRectangle(win_info->display,win_info->screen0,win_info->gc,
			312+2 ,600+80-135+2,
			118-8,135-8);
		XSetForeground(win_info->display,win_info->gc,win_info->foreground); // write point
	}

	// initializ cards on players cards
	for(k=0;k<=4;k++){
		number_of_player_card_x=0;
		number_of_player_card_y=0;
		for(j=1;j<=13;j++){
			for(i=0;i<=3;i++){
				if(players_card[k][i][j]==1){
					xflag=(j-1)*70+2;
					yflag=240+i*80 ;
					XCopyArea(win_info->display, win_info->screen5, win_info->screen0, win_info->gc, 
						xflag,yflag,
						70,80,
						10+70*number_of_player_card_x + 300*seat_x[get_seat(sekijun,k)],
						105+80*number_of_player_card_y +350*seat_y[get_seat(sekijun,k)]);
					number_of_player_card_x++;
					if(number_of_player_card_x>=4){
						number_of_player_card_y++;
						number_of_player_card_x=0;
					}
				}
				if((*accept_flag==k)&&(stage_card[i][j]==1)){
					xflag=(j-1)*70+2;
					yflag=240+i*80 ;
					XCopyArea(win_info->display, win_info->screen5, win_info->screen0, win_info->gc, 
						xflag,yflag,
						70,80,
						10+70*number_of_player_card_x + 300*seat_x[get_seat(sekijun,k)],
						105+80*number_of_player_card_y +350*seat_y[get_seat(sekijun,k)]);
					XSetForeground(win_info->display,win_info->gc,color[players_card[0][5][3]]);
					XSetLineAttributes(win_info->display,win_info->gc,10,LineSolid,CapNotLast,JoinMiter);
					XDrawRectangle(win_info->display,win_info->screen0,win_info->gc,
						10+70*number_of_player_card_x + 300*seat_x[get_seat(sekijun,k)]+4,
						105+80*number_of_player_card_y +350*seat_y[get_seat(sekijun,k)]+4,
						70-8,80-8);
					XSetForeground(win_info->display,win_info->gc,win_info->foreground); // write point
					number_of_player_card_x++;
					if(number_of_player_card_x>=4){
						number_of_player_card_y++;
						number_of_player_card_x=0;
					}
				}
			}
		}
		if(players_card[k][4][1]==2){
			xflag=910+2;
			yflag=240;
			XCopyArea(win_info->display, win_info->screen5, win_info->screen0, win_info->gc, 
				xflag,yflag,
				70,80,
				10+70*number_of_player_card_x + 300*seat_x[get_seat(sekijun,k)],
				105+80*number_of_player_card_y +350*seat_y[get_seat(sekijun,k)]);
	//		XCopyArea(win_info->display, win_info->screen2, win_info->screen0, win_info->gc, xflag,yflag,win_info->card_width,win_info->card_height,win_info->card_offset_x+win_info->card_width*number_of_player_card+win_info->on_stage_player_x_offset,win_info->on_stage_player_size*get_seat(sekijun,k)+win_info->on_stage_player_y_offset);
			number_of_player_card_x++;
		}
		if((*accept_flag==k) && (joker_flag==1)){
			xflag=910+2;
			yflag=240;
			XCopyArea(win_info->display, win_info->screen5, win_info->screen0, win_info->gc, 
				xflag,yflag,
				70,80,
				10+70*number_of_player_card_x + 300*seat_x[get_seat(sekijun,k)],
				105+80*number_of_player_card_y +350*seat_y[get_seat(sekijun,k)]);
			XSetForeground(win_info->display,win_info->gc,color[players_card[0][5][3]]);
			XSetLineAttributes(win_info->display,win_info->gc,10,LineSolid,CapNotLast,JoinMiter);
			XDrawRectangle(win_info->display,win_info->screen0,win_info->gc,
				10+70*number_of_player_card_x + 300*seat_x[get_seat(sekijun,k)]+4,
				105+80*number_of_player_card_y +350*seat_y[get_seat(sekijun,k)]+4,
				70-4,80-4);
			XSetForeground(win_info->display,win_info->gc,win_info->foreground); // write point
		}
	}
	*accept_flag=10;

	// winner mark 
	tmp2=0;
	for(i=0;i<=4;i++){
		if(count_card_num_r(players_card[i],&tmp)==0){
			XCopyArea(win_info->display,win_info->screen6,win_info->screen0,win_info->gc,i*300,0,300,350,300*seat_x[get_seat(sekijun,i)],350*seat_y[get_seat(sekijun,i)]);
			tmp2++;
			XSetForeground(win_info->display,win_info->gc,win_info->background); // clean pixmap
			XFillRectangle(win_info->display, win_info->screen0, win_info->gc,  
				300*seat_x[get_seat(sekijun,i)]+180,350*seat_y[get_seat(sekijun,i)]+300, 
				115, 45);
			XSetForeground(win_info->display,win_info->gc,win_info->foreground);
			XDrawImageString(win_info->display, win_info->screen0, win_info->gc, 
				300*seat_x[get_seat(sekijun,i)]+185,350*seat_y[get_seat(sekijun,i)]+330 , 
				player_name[i], strlen(player_name[i]));
			XCopyArea(win_info->display,win_info->screen5,win_info->screen0,win_info->gc,
				280+100*get_seat(mibun,i),60,
				100,50,
				300*seat_x[get_seat(sekijun,i)]+195,350*seat_y[get_seat(sekijun,i)]+250);
		}else{
			loser=i;
		}
	}
	if(tmp2==4){
		XCopyArea(win_info->display,win_info->screen6,win_info->screen0,win_info->gc,loser*300,350,300,350,300*seat_x[get_seat(sekijun,loser)],350*seat_y[get_seat(sekijun,loser)]);
		XSetForeground(win_info->display,win_info->gc,win_info->background); // clean pixmap
		XFillRectangle(win_info->display, win_info->screen0, win_info->gc,  
			300*seat_x[get_seat(sekijun,loser)]+180,350*seat_y[get_seat(sekijun,loser)]+300, 
			115, 50);
		XSetForeground(win_info->display,win_info->gc,win_info->foreground);
		XDrawImageString(win_info->display, win_info->screen0, win_info->gc, 
			300*seat_x[get_seat(sekijun,loser)]+185,350*seat_y[get_seat(sekijun,loser)]+330 , 
			player_name[loser], strlen(player_name[loser]));
		//XCopyArea(win_info->display,win_info->screen5,win_info->screen0,win_info->gc,
		//	280+100*4,60,
		//	100,50,
		//	300*seat_x[get_seat(sekijun,loser)]+195,350*seat_y[get_seat(sekijun,loser)]+250);
	}

	if(stage_card[5][6]>=100){
		XCopyArea(win_info->display,win_info->screen5,win_info->screen0,win_info->gc,
			0+400* (*accept_flag)%3, 560+200* (*accept_flag>3),
			//400,560,
			400,200,
			250,200);
	}


	XCopyArea(win_info->display, win_info->screen0, win_info->screen3, win_info->gc, 0,0,win_info->window_width, win_info->window_height, 0,0); // delete old card on stage
	XFlush( win_info->display );

}

