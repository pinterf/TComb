.xmm
.model flat,c

.data

align 16

onesByte qword 2 dup(0101010101010101h)
sixsMask_W qword 2 dup(0006000600060006h)
eightsMask_W qword 2 dup(0008000800080008h)

.code

buildFinalMask_SSE2 proc public uses ebx esi edi s1p:dword,s2p:dword,m1p:dword,dstp:dword,stride:dword,width_:dword,height:dword,thresh:dword

	mov eax,s1p
	mov ebx,s2p
	mov edx,m1p
	mov esi,dstp
	mov edi,width_
	
	dec thresh
	movd xmm4,thresh
	punpcklbw xmm4, xmm4
	punpcklwd xmm4, xmm4
	punpckldq xmm4, xmm4
	punpcklqdq xmm4, xmm4
	
	pxor xmm5,xmm5
	
yloop:
	xor ecx,ecx
	align 16
xloop:
	movdqa xmm0,[eax+ecx]
	movdqa xmm1,[ebx+ecx]
	movdqa xmm2,xmm0
	psubusb xmm0,xmm1
	psubusb xmm1,xmm2
	por xmm0,xmm1
	psubusb xmm0,xmm4
	pcmpeqb xmm0,xmm5
	pand xmm0,[edx+ecx]
	movdqa [esi+ecx],xmm0
	
	add ecx,16
	cmp ecx,edi
	jl xloop
	
	add eax,stride
	add ebx,stride
	add edx,stride
	add esi,stride
	dec height
	jnz yloop

	ret

buildFinalMask_SSE2 endp



andNeighborsInPlace_SSE2 proc public uses esi edi srcp:dword,stride:dword,width_:dword,height:dword

	mov eax,srcp
	mov edx,width_
	mov esi,eax
	sub esi,stride
	mov edi,eax
	add edi,stride
	
yloop:
	xor ecx,ecx
	align 16
xloop:
	movdqa xmm0,[esi+ecx]
	movdqu xmm1,[esi+ecx-1]
	por xmm0,xmm1
	movdqu xmm1,[esi+ecx+1]
	por xmm0,xmm1
	movdqa xmm1,[eax+ecx]
	movdqu xmm2,[edi+ecx-1]
	por xmm0,xmm2
	por xmm0,[edi+ecx]
	movdqu xmm2,[edi+ecx+1]
	por xmm0,xmm2
	pand xmm0,xmm1
	movdqa [eax+ecx],xmm0
	
	add ecx,16
	cmp ecx,edx
	jl xloop
	
	add eax,stride
	add esi,stride
	add edi,stride
	dec height
	jnz yloop

	ret

andNeighborsInPlace_SSE2 endp



absDiff_SSE2 proc public uses ebx esi srcp1:dword,srcp2:dword,dstp:dword,stride:dword,width_:dword,height:dword
	
	mov eax,srcp1
	mov esi,srcp2
	mov ebx,dstp
	mov edx,width_
	
yloop:
	xor ecx,ecx
	align 16
xloop:
	movdqa xmm0,[eax+ecx]
	movdqa xmm1,[esi+ecx]
	movdqa xmm2,xmm0
	psubusb xmm0,xmm1
	psubusb xmm1,xmm2
	por xmm0,xmm1
	movdqa [ebx+ecx],xmm0
	
	add ecx,16
	cmp ecx,edx
	jl xloop
	
	add eax,stride
	add esi,stride
	add ebx,stride
	dec height
	jnz yloop

	ret
	
absDiff_SSE2 endp



absDiffAndMinMask_SSE2 proc public uses ebx esi edi srcp1:dword,srcp2:dword,dstp:dword,stride:dword,width_:dword,height:dword
	
	mov eax,srcp1
	mov esi,srcp2
	mov ebx,dstp
	mov edx,width_
	mov edi,height
	
yloop:
	xor ecx,ecx
	align 16
xloop:
	movdqa xmm0,[eax+ecx]
	movdqa xmm1,[esi+ecx]
	movdqa xmm2,xmm0
	psubusb xmm0,xmm1
	psubusb xmm1,xmm2
	por xmm0,xmm1
	pminub xmm0,[ebx+ecx]
	movdqa [ebx+ecx],xmm0
	
	add ecx,16
	cmp ecx,edx
	jl xloop
	
	add eax,stride
	add esi,stride
	add ebx,stride
	dec edi
	jnz yloop

	ret

absDiffAndMinMask_SSE2 endp



absDiffAndMinMaskThresh_SSE2 proc public uses ebx esi edi srcp1:dword,srcp2:dword,dstp:dword,stride:dword,width_:dword,height:dword,thresh:dword
	
	mov eax,srcp1
	mov esi,srcp2
	mov ebx,dstp
	mov edx,width_
	mov edi,height
	
	dec thresh
	movd xmm3,thresh
	punpcklbw xmm3,xmm3
	punpcklwd xmm3,xmm3
	punpckldq xmm3,xmm3
	punpcklqdq xmm3,xmm3
	
	pxor xmm4,xmm4
	
yloop:
	xor ecx,ecx
	align 16
xloop:
	movdqa xmm0,[eax+ecx]
	movdqa xmm1,[esi+ecx]
	movdqa xmm2,xmm0
	psubusb xmm0,xmm1
	psubusb xmm1,xmm2
	por xmm0,xmm1
	pminub xmm0,[ebx+ecx]
	psubusb xmm0,xmm3
	pcmpeqb xmm0,xmm4
	movdqa [ebx+ecx],xmm0
	
	add ecx,16
	cmp ecx,edx
	jl xloop
	
	add eax,stride
	add esi,stride
	add ebx,stride
	dec edi
	jnz yloop

	ret
	
absDiffAndMinMaskThresh_SSE2 endp



MinMax_SSE2 proc public uses ebx esi edi srcp:dword,minp:dword,maxp:dword,src_stride:dword,min_stride:dword,width_:dword,height:dword,thresh:dword
	
	mov eax,srcp
	mov esi,eax
	sub esi,src_stride
	mov edi,eax
	add edi,src_stride
	mov ebx,minp
	mov edx,maxp
	
	movd xmm3,thresh
	punpcklbw xmm3,xmm3
	punpcklwd xmm3,xmm3
	punpckldq xmm3,xmm3
	punpcklqdq xmm3,xmm3
	
yloop:
	xor ecx,ecx
	align 16
xloop:
	; srcp-1 is aligned because the pointer passed to this function is srcp+stride+1.
	movdqa xmm0,[esi+ecx-1]
	movdqa xmm1,xmm0
	movdqu xmm2,[esi+ecx]
	pminub xmm0,xmm2
	pmaxub xmm1,xmm2
	movdqu xmm2,[esi+ecx+1]
	pminub xmm0,xmm2
	pmaxub xmm1,xmm2
	movdqa xmm2,[eax+ecx-1]
	pminub xmm0,xmm2
	pmaxub xmm1,xmm2
	movdqu xmm2,[eax+ecx]
	pminub xmm0,xmm2
	pmaxub xmm1,xmm2
	movdqu xmm2,[eax+ecx+1]
	pminub xmm0,xmm2
	pmaxub xmm1,xmm2
	movdqa xmm2,[edi+ecx-1]
	pminub xmm0,xmm2
	pmaxub xmm1,xmm2
	movdqu xmm2,[edi+ecx]
	pminub xmm0,xmm2
	pmaxub xmm1,xmm2
	movdqu xmm2,[edi+ecx+1]
	pminub xmm0,xmm2
	pmaxub xmm1,xmm2
	psubusb xmm0,xmm3
	paddusb xmm1,xmm3
	movdqa [ebx+ecx],xmm0
	movdqa [edx+ecx],xmm1
	
	add ecx,16
	cmp ecx,width_
	jl xloop
	
	add esi,src_stride
	add eax,src_stride
	add edi,src_stride
	add ebx,min_stride
	add edx,min_stride
	dec height
	jnz yloop

	ret
	
MinMax_SSE2 endp



checkOscillation5_SSE2 proc public uses ebx esi edi p2p:dword,p1p:dword,s1p:dword,n1p:dword,n2p:dword,dstp:dword,stride:dword,width_:dword,height:dword,thresh:dword

	mov eax,p2p
	mov ebx,p1p
	mov edx,s1p
	mov edi,n1p
	mov esi,n2p
	
	
	pxor xmm6,xmm6
	
	dec thresh
	movd xmm7,thresh
	punpcklbw xmm7,xmm7
	punpcklwd xmm7,xmm7
	punpckldq xmm7,xmm7
	punpcklqdq xmm7,xmm7
	
yloop:
	xor ecx,ecx
	align 16
xloop:
	movdqa xmm0,[eax+ecx]
	movdqa xmm2,[ebx+ecx]
	movdqa xmm1,xmm0
	movdqa xmm3,xmm2
	pminub xmm0,[edx+ecx]
	pmaxub xmm1,[edx+ecx]
	pminub xmm2,[edi+ecx]
	pmaxub xmm3,[edi+ecx]
	pminub xmm0,[esi+ecx]
	pmaxub xmm1,[esi+ecx]
	movdqa xmm4,xmm3
	movdqa xmm5,xmm1
	psubusb xmm4,xmm2
	psubusb xmm5,xmm0
	psubusb xmm4,xmm7
	psubusb xmm5,xmm7
	psubusb xmm2,oword ptr onesByte
	psubusb xmm0,oword ptr onesByte
	psubusb xmm1,xmm2
	psubusb xmm3,xmm0
	pcmpeqb xmm1,xmm6
	pcmpeqb xmm3,xmm6
	pcmpeqb xmm4,xmm6
	pcmpeqb xmm5,xmm6
	mov eax,dstp
	por xmm1,xmm3
	pand xmm4,xmm5
	pand xmm1,xmm4
	movdqa [eax+ecx],xmm1
	
	add ecx,16
	mov eax,p2p
	cmp ecx,width_
	jl xloop
	
	mov eax,stride
	add ebx,stride
	add p2p,eax
	add edx,stride
	add edi,stride
	add dstp,eax
	add esi,stride
	mov eax,p2p
	dec height
	jnz yloop

	ret
	
checkOscillation5_SSE2 endp



calcAverages_SSE2 proc public uses ebx esi edi s1p:dword,s2p:dword,dstp:dword,stride:dword,width_:dword,height:dword
	
	mov eax,s1p
	mov ebx,s2p
	mov edx,dstp
	mov edi,height
	mov esi,width_
	
yloop:
	xor ecx,ecx
	align 16
xloop:
	movdqa xmm0,[eax+ecx]
	pavgb xmm0,[ebx+ecx]
	movdqa [edx+ecx],xmm0
	
	add ecx,16
	cmp ecx,esi
	jl xloop
	
	add eax,stride
	add ebx,stride
	add edx,stride
	dec edi
	jnz yloop

	ret
	
calcAverages_SSE2 endp



checkAvgOscCorrelation_SSE2 proc public uses ebx esi edi s1p:dword,s2p:dword,s3p:dword,s4p:dword,dstp:dword,stride:dword,width_:dword,height:dword,thresh:dword
	
	mov eax,s1p
	mov ebx,s2p
	mov edx,s3p
	mov edi,s4p
	mov esi,dstp
	
	dec thresh
	movd xmm2, thresh
	punpcklbw xmm2, xmm2
	punpcklwd xmm2, xmm2
	punpckldq xmm2, xmm2
	punpcklqdq xmm2, xmm2
	
	pxor xmm3,xmm3
	
yloop:
	xor ecx,ecx
	align 16
xloop:
	movdqa xmm5,[eax+ecx]
	movdqa xmm0,xmm5
	movdqa xmm1,xmm5
	movdqa xmm5,[ebx+ecx]
	pminub xmm0,xmm5
	pmaxub xmm1,xmm5
	movdqa xmm5,[edx+ecx]
	pminub xmm0,xmm5
	pmaxub xmm1,xmm5
	movdqa xmm5,[edi+ecx]
	pminub xmm0,xmm5
	pmaxub xmm1,xmm5
	psubusb xmm1,xmm0
	movdqa xmm4,[esi+ecx]
	psubusb xmm1,xmm2
	pcmpeqb xmm1,xmm3
	pand xmm1,xmm4
	movdqa [esi+ecx],xmm1
	
	add ecx,16
	cmp ecx,width_
	jl xloop
	
	add eax,stride
	add ebx,stride
	add edx,stride
	add edi,stride
	add esi,stride
	dec height
	jnz yloop

	ret
	
checkAvgOscCorrelation_SSE2 endp



or3Masks_SSE2 proc public uses ebx esi edi s1p:dword,s2p:dword,s3p:dword,dstp:dword,stride:dword,width_:dword,height:dword
	
	mov eax,s1p
	mov ebx,s2p
	mov edx,s3p
	mov edi,dstp
	mov esi,width_
	
yloop:
	xor ecx,ecx
	align 16
xloop:
	movdqa xmm0,[eax+ecx]
	por xmm0,[ebx+ecx]
	por xmm0,[edx+ecx]
	movdqa [edi+ecx],xmm0
	
	add ecx,16
	cmp ecx,esi
	jl xloop
	
	add eax,stride
	add ebx,stride
	add edx,stride
	add edi,stride
	dec height
	jnz yloop

	ret
	
or3Masks_SSE2 endp



orAndMasks_SSE2 proc public uses ebx esi edi s1p:dword,s2p:dword,dstp:dword,stride:dword,width_:dword,height:dword
	
	mov eax,s1p
	mov ebx,s2p
	mov edx,dstp
	mov edi,width_
	mov esi,height
	
yloop:
	xor ecx,ecx
	align 16
xloop:
	movdqa xmm0,[eax+ecx]
	movdqa xmm1,[edx+ecx]
	pand xmm0,[ebx+ecx]
	por xmm1,xmm0
	movdqa [edx+ecx],xmm1
	
	add ecx,16
	cmp ecx,edi
	jl xloop
	
	add eax,stride
	add ebx,stride
	add edx,stride
	dec esi
	jnz yloop

	ret
	
orAndMasks_SSE2 endp



andMasks_SSE2 proc public uses ebx esi edi s1p:dword,s2p:dword,dstp:dword,stride:dword,width_:dword,height:dword
	
	mov eax,s1p
	mov ebx,s2p
	mov edx,dstp
	mov edi,width_
	mov esi,height
	
yloop:
	xor ecx,ecx
	align 16
xloop:
	movdqa xmm0,[eax+ecx]
	pand xmm0,[ebx+ecx]
	movdqa [edx+ecx],xmm0
	
	add ecx,16
	cmp ecx,edi
	jl xloop
	
	add eax,stride
	add ebx,stride
	add edx,stride
	dec esi
	jnz yloop

	ret
	
andMasks_SSE2 endp



checkSceneChange_SSE2 proc public uses ebx esi edi s1p:dword,s2p:dword,stride:dword,width_:dword,height:dword,diffp:dword
	
	mov eax,s1p
	mov edi,s2p
	mov esi,stride
	mov edx,width_
	pxor xmm1,xmm1
	
yloop:
	xor ecx,ecx
	align 16
xloop:
	movdqa xmm0,[eax+ecx]
	psadbw xmm0,[edi+ecx]
	paddd xmm1,xmm0

	add ecx,16
	cmp ecx,edx
	jl xloop

	add eax,esi
	add edi,esi
	dec height
	jnz yloop

	movdqa xmm2,xmm1
	psrldq xmm1,8
	paddd xmm2,xmm1

	mov eax, diffp
	movd DWORD PTR [eax],xmm2

	ret
	
checkSceneChange_SSE2 endp



VerticalBlur3_SSE2 proc public uses ebx esi edi srcp:dword,dstp:dword,stride:dword,width_:dword,height:dword
	
	mov eax,srcp
	mov ebx,dstp
	mov edx,stride
	mov esi,eax
	mov edi,eax
	sub esi,edx
	add edi,edx
	mov edx,width_
	
	; 0x0002,for rounding
	pcmpeqb xmm6,xmm6
	psrlw xmm6,15
	psllw xmm6,1
	pxor xmm7,xmm7

	xor ecx,ecx

toploop:
	movdqa xmm0,[eax+ecx]
	pavgb xmm0,[edi+ecx]
	movdqa [ebx+ecx],xmm0
	
	add ecx,16
	cmp ecx,edx
	jl toploop
	
	add esi,stride
	add eax,stride
	add edi,stride
	add ebx,stride
	sub height,2 ; the main loop processes 2 lines fewer than the height
	
yloop:
	xor ecx,ecx
xloop:
	movdqa xmm0,[esi+ecx]
	movdqa xmm1,[eax+ecx]
	movdqa xmm2,[edi+ecx]
	movdqa xmm3,xmm0
	movdqa xmm4,xmm1
	movdqa xmm5,xmm2
	punpcklbw xmm0,xmm7
	punpcklbw xmm1,xmm7
	punpcklbw xmm2,xmm7
	punpckhbw xmm3,xmm7
	punpckhbw xmm4,xmm7
	punpckhbw xmm5,xmm7
	
	; add bottom to top
	paddw xmm0,xmm2
	paddw xmm3,xmm5
	
	; multiply center by 2
	psllw xmm1,1
	psllw xmm4,1
	
	; add center to sum
	paddw xmm0,xmm1
	paddw xmm3,xmm4
	
	; add 2 to sum
	paddw xmm0,xmm6
	paddw xmm3,xmm6
	
	; divide by 4
	psrlw xmm0,2
	psrlw xmm3,2
	packuswb xmm0,xmm3
	movdqa [ebx+ecx],xmm0
	
	add ecx,16
	cmp ecx,edx
	jl xloop
	
	add esi,stride
	add eax,stride
	add edi,stride
	add ebx,stride
	dec height
	jnz yloop
	
	xor ecx,ecx
	
bottomloop:
	movdqa xmm0,[esi+ecx]
	pavgb xmm0,[eax+ecx]
	movdqa [ebx+ecx],xmm0
	
	add ecx,16
	cmp ecx,edx
	jl bottomloop

	ret
	
VerticalBlur3_SSE2 endp



HorizontalBlur3_SSE2 proc public srcp:dword,dstp:dword,stride:dword,width_:dword,height:dword
	
	mov eax,srcp
	mov edx,dstp
	pxor xmm7,xmm7
	; 0x0002,for rounding
	pcmpeqb xmm6,xmm6
	psrlw xmm6,15
	psllw xmm6,1

yloop:
	xor ecx,ecx
	align 16
xloop:
	movdqu xmm0,[eax+ecx-1]
	movdqa xmm1,[eax+ecx]
	movdqu xmm2,[eax+ecx+1]
	movdqa xmm3,xmm0
	movdqa xmm4,xmm1
	movdqa xmm5,xmm2
	punpcklbw xmm0,xmm7
	punpcklbw xmm1,xmm7
	punpcklbw xmm2,xmm7
	punpckhbw xmm3,xmm7
	punpckhbw xmm4,xmm7
	punpckhbw xmm5,xmm7
	; center * 2
	psllw xmm1,1
	psllw xmm4,1
	paddw xmm1,xmm0
	paddw xmm4,xmm3
	paddw xmm1,xmm2
	paddw xmm4,xmm5

	; add 2 to sum
	paddw xmm1,xmm6
	paddw xmm4,xmm6

	; divide by 4
	psrlw xmm1,2
	psrlw xmm4,2
	packuswb xmm1,xmm4
	movdqa [edx+ecx],xmm1
	
	add ecx,16
	cmp ecx,width_
	jl xloop
	
	add eax,stride
	add edx,stride
	dec height
	jnz yloop

	ret
	
HorizontalBlur3_SSE2 endp



HorizontalBlur6_SSE2 proc public srcp:dword,dstp:dword,stride:dword,width_:dword,height:dword
	
	mov eax,srcp
	mov edx,dstp
	movdqu xmm6,oword ptr sixsMask_W
	pxor xmm7,xmm7
	
yloop:
	xor ecx,ecx
	align 16
xloop:
	movdqu xmm0,[eax+ecx-2]
	movdqu xmm1,[eax+ecx+2]
	movdqa xmm2,xmm0
	movdqa xmm3,xmm1
	punpcklbw xmm0,xmm7
	punpcklbw xmm1,xmm7
	punpckhbw xmm2,xmm7
	punpckhbw xmm3,xmm7
	
	; srcp[x-2] + srcp[x+2]
	paddw xmm0,xmm1
	paddw xmm2,xmm3
	
	; srcp[x-1] + srcp[x+1]
	movdqu xmm1,[eax+ecx-1]
	movdqu xmm3,[eax+ecx+1]
	movdqa xmm4,xmm1
	movdqa xmm5,xmm3
	punpcklbw xmm1,xmm7
	punpcklbw xmm3,xmm7
	punpckhbw xmm4,xmm7
	punpckhbw xmm5,xmm7
	paddw xmm1,xmm3
	paddw xmm4,xmm5
	
	; (srcp[x-1 + srcp[x+])*4
	psllw xmm1,2
	psllw xmm4,2
	
	; (srcp[x-1 + srcp[x+])*4 + srcp[x-2] + srcp[x+2]
	paddw xmm0,xmm1
	paddw xmm2,xmm4
	
	; srcp[x] * 6
	movdqa xmm1,[eax+ecx]
	movdqu xmm5,oword ptr eightsMask_W
	movdqa xmm3,xmm1
	punpcklbw xmm1,xmm7
	punpckhbw xmm3,xmm7
	pmullw xmm1,xmm6
	pmullw xmm3,xmm6
	paddw xmm0,xmm1
	paddw xmm2,xmm3
	
	; add 8
	paddw xmm0,xmm5
	paddw xmm2,xmm5
	
	; divide by 16
	psrlw xmm0,4
	psrlw xmm2,4
	packuswb xmm0,xmm2
	movdqa [edx+ecx],xmm0
	
	add ecx,16
	cmp ecx,width_
	jl xloop
	
	add eax,stride
	add edx,stride
	dec height
	jnz yloop

	ret
	
HorizontalBlur6_SSE2 endp



end