.code

;buildFinalMask_SSE2 proc s1p:dword,s2p:dword,m1p:dword,dstp:dword,stride:dword,width_:dword,height:dword,thresh:dword
; s1p = rcx
; s2p = rdx
; m1p = r8
; dstp = r9

buildFinalMask_SSE2 proc public frame
	
	stride equ dword ptr [rbp+48]
	width_ equ dword ptr [rbp+56]
	height equ dword ptr [rbp+64]
	thresh equ dword ptr [rbp+72]
	
	push rbp
.pushreg rbp
	mov rbp,rsp
	push rbx
.pushreg rbx
	push rsi
.pushreg rsi
	push rdi
.pushreg rdi
.endprolog
	
	mov rax,rcx
	mov rbx,rdx
	mov rdx,r8
	mov rsi,r9
	movsxd r8,stride
	xor rdi,rdi
	mov edi,width_
	xor r9,r9
	mov r9d,height
	mov r10,16

	dec thresh
	movd xmm4,thresh
	punpcklbw xmm4, xmm4
	punpcklwd xmm4, xmm4
	punpckldq xmm4, xmm4
	punpcklqdq xmm4, xmm4
	pxor xmm5,xmm5

yloop:
	xor rcx,rcx
xloop:
	movdqa xmm0,[rax+rcx]
	movdqa xmm1,[rbx+rcx]
	movdqa xmm2,xmm0
	psubusb xmm0,xmm1
	psubusb xmm1,xmm2
	por xmm0,xmm1
	psubusb xmm0,xmm4
	pcmpeqb xmm0,xmm5
	pand xmm0,[rdx+rcx]
	movdqa [rsi+rcx],xmm0

	add rcx,r10
	cmp rcx,rdi
	jl xloop

	add rax,r8
	add rbx,r8
	add rdx,r8
	add rsi,r8
	dec r9
	jnz yloop
	
	pop rdi
	pop rsi
	pop rbx
	pop rbp

	ret

buildFinalMask_SSE2 endp



;andNeighborsInPlace_SSE2 proc srcp:dword,stride:dword,width_:dword,height:dword
; srcp = rcx
; stride = rdx
; width_ = r8d
; height = r9d

andNeighborsInPlace_SSE2 proc public frame
	
	push rbp
.pushreg rbp
	mov rbp,rsp
	push rsi
.pushreg rsi
	push rdi
.pushreg rdi
.endprolog

	mov rax,rcx
	xchg r8,rdx
	movsxd r8,r8d
	mov rsi,rax
	mov rdi,rax
	sub rsi,r8
	add rdi,r8
	mov r10,16

yloop:
	xor rcx,rcx
xloop:
	movdqa xmm0,[rsi+rcx]
	movdqu xmm1,[rsi+rcx-1]
	por xmm0,xmm1
	movdqu xmm1,[rsi+rcx+1]
	por xmm0,xmm1
	movdqa xmm1,[rax+rcx]
	movdqu xmm2,[rdi+rcx-1]
	por xmm0,xmm2
	por xmm0,[rdi+rcx]
	movdqu xmm2,[rdi+rcx+1]
	por xmm0,xmm2
	pand xmm0,xmm1
	movdqa [rax+rcx],xmm0

	add rcx,r10
	cmp rcx,rdx
	jl xloop

	add rax,r8
	add rsi,r8
	add rdi,r8
	dec r9d
	jnz yloop

	pop rdi
	pop rsi
	pop rbp
	
	ret

andNeighborsInPlace_SSE2 endp



;absDiff_SSE2 proc srcp1:dword,srcp2:dword,dstp:dword,stride:dword,width_:dword,height:dword
; srcp1 = rcx
; srcp2 = rdx
; dstp = r8
; stride = r9d

absDiff_SSE2 proc public frame

	width_ equ dword ptr [rbp+48]
	height equ dword ptr [rbp+56]
	
	push rbp
.pushreg rbp
	mov rbp,rsp
	push rbx
.pushreg rbx
	push rsi
.pushreg rsi
.endprolog
	
	mov rax,rcx
	mov rsi,rdx
	mov rbx,r8
	movsxd r8,r9d
	xor rdx,rdx
	mov edx,width_
	xor r9,r9
	mov r9d,height
	mov r10,16
	
yloop:
	xor rcx,rcx
xloop:
	movdqa xmm0,[rax+rcx]
	movdqa xmm1,[rsi+rcx]
	movdqa xmm2,xmm0
	psubusb xmm0,xmm1
	psubusb xmm1,xmm2
	por xmm0,xmm1
	movdqa [rbx+rcx],xmm0
	
	add rcx,r10
	cmp rcx,rdx
	jl xloop
	
	add rax,r8
	add rsi,r8
	add rbx,r8
	dec r9d
	jnz yloop
	
	pop rsi
	pop rbx
	pop rbp

	ret
	
absDiff_SSE2 endp



;absDiffAndMinMask_SSE2 proc srcp1:dword,srcp2:dword,dstp:dword,stride:dword,width_:dword,height:dword
; srcp1 = rcx
; srcp2 = rdx
; dstp = r8
; stride = r9d

absDiffAndMinMask_SSE2 proc public frame

	width_ equ dword ptr [rbp+48]
	height equ dword ptr [rbp+56]
	
	push rbp
.pushreg rbp
	mov rbp,rsp
	push rbx
.pushreg rbx
	push rsi
.pushreg rsi
	push rdi
.pushreg rdi
.endprolog
	
	mov rax,rcx
	mov rsi,rdx
	mov rbx,r8
	movsxd r8,r9d
	xor rdx,rdx
	mov edx,width_
	xor rdi,rdi
	mov edi,height
	mov r10,16
	
yloop:
	xor rcx,rcx
xloop:
	movdqa xmm0,[rax+rcx]
	movdqa xmm1,[rsi+rcx]
	movdqa xmm2,xmm0
	psubusb xmm0,xmm1
	psubusb xmm1,xmm2
	por xmm0,xmm1
	pminub xmm0,[rbx+rcx]
	movdqa [rbx+rcx],xmm0
	
	add rcx,r10
	cmp rcx,rdx
	jl xloop
	
	add rax,r8
	add rsi,r8
	add rbx,r8
	dec edi
	jnz yloop

	pop rdi
	pop rsi
	pop rbx
	pop rbp
	
	ret

absDiffAndMinMask_SSE2 endp



;absDiffAndMinMaskThresh_SSE2 proc srcp1:dword,srcp2:dword,dstp:dword,stride:dword,width_:dword,height:dword,thresh:dword
; srcp1 = rcx
; srcp2 = rdx
; dstp = r8
; stride = r9d

absDiffAndMinMaskThresh_SSE2 proc public frame

	width_ equ dword ptr [rbp+48]
	height equ dword ptr [rbp+56]
	thresh equ dword ptr [rbp+64]
	
	push rbp
.pushreg rbp
	mov rbp,rsp
	push rbx
.pushreg rbx
	push rsi
.pushreg rsi
	push rdi
.pushreg rdi
.endprolog
	
	mov rax,rcx
	mov rsi,rdx
	mov rbx,r8
	movsxd r8,r9d
	xor rdx,rdx
	mov edx,width_
	xor rdi,rdi
	mov edi,height
	dec thresh
	movd xmm3,thresh
	punpcklbw xmm3,xmm3
	punpcklwd xmm3,xmm3
	punpckldq xmm3,xmm3
	punpcklqdq xmm3,xmm3
	pxor xmm4,xmm4
	mov r10,16
	
yloop:
	xor rcx,rcx
xloop:
	movdqa xmm0,[rax+rcx]
	movdqa xmm1,[rsi+rcx]
	movdqa xmm2,xmm0
	psubusb xmm0,xmm1
	psubusb xmm1,xmm2
	por xmm0,xmm1
	pminub xmm0,[rbx+rcx]
	psubusb xmm0,xmm3
	pcmpeqb xmm0,xmm4
	movdqa [rbx+rcx],xmm0
	
	add rcx,r10
	cmp rcx,rdx
	jl xloop
	
	add rax,r8
	add rsi,r8
	add rbx,r8
	dec edi
	jnz yloop
	
	pop rdi
	pop rsi
	pop rbx
	pop rbp

	ret
	
absDiffAndMinMaskThresh_SSE2 endp



;MinMax_SSE2 proc srcp:dword,minp:dword,maxp:dword,src_stride:dword,min_stride:dword,width_:dword,height:dword,thresh:dword
; srcp = rcx
; minp = edx
; maxp = r8d
; src_stride = r9d

MinMax_SSE2 proc public frame

	min_stride equ dword ptr [rbp+48]
	width_ equ dword ptr [rbp+56]
	height equ dword ptr [rbp+64]
	thresh equ dword ptr [rbp+72]
	
	push rbp
.pushreg rbp
	mov rbp,rsp
	push rbx
.pushreg rbx
	push rsi
.pushreg rsi
	push rdi
.pushreg rdi
	push r12
.pushreg r12
.endprolog
	
	mov rax,rcx
	mov rsi,rax
	mov rdi,rax
	mov rbx,rdx
	mov rdx,r8
	movsxd r8,r9d
	movsxd r9,min_stride
	mov r10d,width_
	mov r11d,height	
	mov r12,16
	sub rsi,r8
	add rdi,r8
	
	movd xmm3,thresh
	punpcklbw xmm3,xmm3
	punpcklwd xmm3,xmm3
	punpckldq xmm3,xmm3
	punpcklqdq xmm3,xmm3
	
yloop:
	xor rcx,rcx
xloop:
	; srcp-1 is aligned because the pointer passed to this function is srcp+stride+1.
	movdqa xmm0,[rsi+rcx-1]
	movdqa xmm1,xmm0
	movdqu xmm2,[rsi+rcx]
	pminub xmm0,xmm2
	pmaxub xmm1,xmm2
	movdqu xmm2,[rsi+rcx+1]
	pminub xmm0,xmm2
	pmaxub xmm1,xmm2
	movdqa xmm2,[rax+rcx-1]
	pminub xmm0,xmm2
	pmaxub xmm1,xmm2
	movdqu xmm2,[rax+rcx]
	pminub xmm0,xmm2
	pmaxub xmm1,xmm2
	movdqu xmm2,[rax+rcx+1]
	pminub xmm0,xmm2
	pmaxub xmm1,xmm2
	movdqa xmm2,[rdi+rcx-1]
	pminub xmm0,xmm2
	pmaxub xmm1,xmm2
	movdqu xmm2,[rdi+rcx]
	pminub xmm0,xmm2
	pmaxub xmm1,xmm2
	movdqu xmm2,[rdi+rcx+1]
	pminub xmm0,xmm2
	pmaxub xmm1,xmm2
	psubusb xmm0,xmm3
	paddusb xmm1,xmm3
	movdqa [rbx+rcx],xmm0
	movdqa [rdx+rcx],xmm1
	
	add rcx,r12
	cmp rcx,r10
	jl xloop
	
	add rsi,r8
	add rax,r8
	add rdi,r8
	add rbx,r9
	add rdx,r9	
	dec r11d
	jnz yloop
	
	pop r12
	pop rdi
	pop rsi
	pop rbx
	pop rbp

	ret
	
MinMax_SSE2 endp



;checkOscillation5_SSE2 proc p2p:dword,p1p:dword,s1p:dword,n1p:dword,n2p:dword,dstp:dword,stride:dword,width_:dword,height:dword,thresh:dword
; p2p = rcx
; p1p = rdx
; s1p = r8
; n1p = r9

checkOscillation5_SSE2 proc public frame

	n2p equ qword ptr [rbp+48]
	dstp equ qword ptr [rbp+56]
	stride equ dword ptr [rbp+64]
	width_ equ dword ptr [rbp+72]
	height equ dword ptr [rbp+80]
	thresh equ dword ptr [rbp+88]
	
	push rbp
.pushreg rbp
	mov rbp,rsp
	push rbx
.pushreg rbx
	push rsi
.pushreg rsi
	push rdi
.pushreg rdi
	sub rsp,64
.allocstack 64
	movdqu oword ptr[rsp],xmm6
.savexmm128 xmm6,0
	movdqu oword ptr[rsp+16],xmm7
.savexmm128 xmm7,16
	movdqu oword ptr[rsp+32],xmm8
.savexmm128 xmm8,32
	movdqu oword ptr[rsp+48],xmm9
.savexmm128 xmm9,48
.endprolog
	
	mov rax,rcx
	mov rbx,rdx
	mov rdx,r8
	mov rdi,r9
	mov rsi,n2p
	mov r8,dstp
	movsxd r9,stride
	mov r10d,width_
	mov r11d,height
	mov r12,16
	
	pxor xmm6,xmm6
	
	dec thresh
	movd xmm7,thresh
	punpcklbw xmm7,xmm7
	punpcklwd xmm7,xmm7
	punpckldq xmm7,xmm7
	punpcklqdq xmm7,xmm7
	
	pcmpeqb xmm9,xmm9
	psrlw xmm9,15
	movdqa xmm8,xmm9
	psllw xmm8,8
	por xmm9,xmm8

yloop:
	xor rcx,rcx
xloop:
	movdqa xmm0,[rax+rcx]
	movdqa xmm2,[rbx+rcx]
	movdqa xmm1,xmm0
	movdqa xmm3,xmm2
	movdqa xmm8,[rdx+rcx]
	pminub xmm0,xmm8
	pmaxub xmm1,xmm8
	movdqa xmm8,[rdi+rcx]
	pminub xmm2,xmm8
	pmaxub xmm3,xmm8
	movdqa xmm8,[rsi+rcx]
	pminub xmm0,xmm8
	pmaxub xmm1,xmm8
	movdqa xmm4,xmm3
	movdqa xmm5,xmm1
	psubusb xmm4,xmm2
	psubusb xmm5,xmm0
	psubusb xmm4,xmm7
	psubusb xmm5,xmm7
	psubusb xmm2,xmm9
	psubusb xmm0,xmm9
	psubusb xmm1,xmm2
	psubusb xmm3,xmm0
	pcmpeqb xmm1,xmm6
	pcmpeqb xmm3,xmm6
	pcmpeqb xmm4,xmm6
	pcmpeqb xmm5,xmm6
	por xmm1,xmm3
	pand xmm4,xmm5
	pand xmm1,xmm4
	movdqa [r8+rcx],xmm1

	add rcx,r12
	cmp rcx,r10
	jl xloop

	add rax,r9
	add rbx,r9
	add rdx,r9
	add rdi,r9
	add rsi,r9
	add r8,r9
	dec r11d
	jnz yloop
	
	movdqu xmm9,oword ptr[rsp+48]
	movdqu xmm8,oword ptr[rsp+32]
	movdqu xmm7,oword ptr[rsp+16]
	movdqu xmm6,oword ptr[rsp]
	add rsp,64
	pop rdi
	pop rsi
	pop rbx
	pop rbp

	ret
	
checkOscillation5_SSE2 endp



;calcAverages_SSE2 proc s1p:dword,s2p:dword,dstp:dword,stride:dword,width_:dword,height:dword
; s1p = rcx
; s2p = rdx
; dstp = r8
; stride = r9d

calcAverages_SSE2 proc public frame
	
	width_ equ dword ptr [rbp+48]
	height equ dword ptr [rbp+56]

	push rbp
.pushreg rbp
	mov rbp,rsp
	push rbx
.pushreg rbx
	push rsi
.pushreg rsi
	push rdi
.pushreg rdi
.endprolog

	mov rax,rcx
	mov rbx,rdx
	mov rdx,r8
	movsxd r8,r9d
	xor rdi,rdi
	mov edi,height
	xor rsi,rsi
	mov esi,width_
	mov r10,16

yloop:
	xor rcx,rcx
xloop:
	movdqa xmm0,[rax+rcx]
	pavgb xmm0,[rbx+rcx]
	movdqa [rdx+rcx],xmm0

	add rcx,r10
	cmp rcx,rsi
	jl xloop

	add rax,r8
	add rbx,r8
	add rdx,r8
	dec edi
	jnz yloop

	pop rdi
	pop rsi
	pop rbx
	pop rbp

	ret
	
calcAverages_SSE2 endp



;checkAvgOscCorrelation_SSE2 proc s1p:dword,s2p:dword,s3p:dword,s4p:dword,dstp:dword,stride:dword,width_:dword,height:dword,thresh:dword
; s1p = rcx
; s2p = rdx
; s3p = r8
; s4p = r9

checkAvgOscCorrelation_SSE2 proc public frame

	dstp equ qword ptr [rbp+48]
	stride equ dword ptr [rbp+56]
	width_ equ dword ptr [rbp+64]
	height equ dword ptr [rbp+72]
	thresh equ dword ptr [rbp+80]

	push rbp
.pushreg rbp
	mov rbp,rsp
	push rbx
.pushreg rbx
	push rsi
.pushreg rsi
	push rdi
.pushreg rdi
.endprolog

	mov rax,rcx
	mov rbx,rdx
	mov rdx,r8
	mov rdi,r9
	mov rsi,dstp
	movsxd r8,stride
	xor r9,r9
	mov r9d,width_
	mov r10d,height
	mov r11,16

	dec thresh
	movd xmm2, thresh
	punpcklbw xmm2, xmm2
	punpcklwd xmm2, xmm2
	punpckldq xmm2, xmm2
	punpcklqdq xmm2, xmm2
	
	pxor xmm3,xmm3

yloop:
	xor rcx,rcx
xloop:
	movdqa xmm5,[rax+rcx]
	movdqa xmm0,xmm5
	movdqa xmm1,xmm5
	movdqa xmm5,[rbx+rcx]
	pminub xmm0,xmm5
	pmaxub xmm1,xmm5
	movdqa xmm5,[rdx+rcx]
	pminub xmm0,xmm5
	pmaxub xmm1,xmm5
	movdqa xmm5,[rdi+rcx]
	pminub xmm0,xmm5
	pmaxub xmm1,xmm5
	psubusb xmm1,xmm0
	movdqa xmm4,[rsi+rcx]
	psubusb xmm1,xmm2
	pcmpeqb xmm1,xmm3
	pand xmm1,xmm4
	movdqa [rsi+rcx],xmm1

	add rcx,r11
	cmp rcx,r9
	jl xloop

	add rax,r8
	add rbx,r8
	add rdx,r8
	add rdi,r8
	add rsi,r8
	dec r10d
	jnz yloop

	pop rdi
	pop rsi
	pop rbx
	pop rbp

	ret
	
checkAvgOscCorrelation_SSE2 endp



;or3Masks_SSE2 proc s1p:dword,s2p:dword,s3p:dword,dstp:dword,stride:dword,width_:dword,height:dword
; s1p = rcx
; s2p = rdx
; s3p = r8
; dstp = r9

or3Masks_SSE2 proc public frame
	
	stride equ dword ptr [rbp+48]
	width_ equ dword ptr [rbp+56]
	height equ dword ptr [rbp+64]

	push rbp
.pushreg rbp
	mov rbp,rsp
	push rbx
.pushreg rbx
	push rsi
.pushreg rsi
	push rdi
.pushreg rdi
.endprolog

	mov rax,rcx
	mov rbx,rdx
	mov rdx,r8
	mov rdi,r9
	movsxd r8,stride
	xor rsi,rsi
	mov esi,width_
	xor r9,r9
	mov r9d,height
	mov r10,16

yloop:
	xor rcx,rcx
xloop:
	movdqa xmm0,[rax+rcx]
	por xmm0,[rbx+rcx]
	por xmm0,[rdx+rcx]
	movdqa [rdi+rcx],xmm0

	add rcx,r10
	cmp rcx,rsi
	jl xloop

	add rax,r8
	add rbx,r8
	add rdx,r8
	add rdi,r8
	dec r9d
	jnz yloop

	pop rdi
	pop rsi
	pop rbx
	pop rbp

	ret
	
or3Masks_SSE2 endp



;orAndMasks_SSE2 proc s1p:dword,s2p:dword,dstp:dword,stride:dword,width_:dword,height:dword
; s1p = rcx
; s2p = rdx
; dstp = r8
; stride = r9d

orAndMasks_SSE2 proc public frame

	width_ equ dword ptr [rbp+48]
	height equ dword ptr [rbp+56]

	push rbp
.pushreg rbp
	mov rbp,rsp
	push rbx
.pushreg rbx
	push rsi
.pushreg rsi
	push rdi
.pushreg rdi
.endprolog

	mov rax,rcx
	mov rbx,rdx
	mov rdx,r8
	movsxd r8,r9d
	xor rdi,rdi
	mov edi,width_
	xor rsi,rsi
	mov esi,height
	mov r10,16

yloop:
	xor rcx,rcx
xloop:
	movdqa xmm0,[rax+rcx]
	movdqa xmm1,[rdx+rcx]
	pand xmm0,[rbx+rcx]
	por xmm1,xmm0
	movdqa [rdx+rcx],xmm1
	add rcx,16
	cmp rcx,rdi
	jl xloop

	add rax,r8
	add rbx,r8
	add rdx,r8
	dec esi
	jnz yloop

	pop rdi
	pop rsi
	pop rbx
	pop rbp

	ret
	
orAndMasks_SSE2 endp



;andMasks_SSE2 proc s1p:dword,s2p:dword,dstp:dword,stride:dword,width_:dword,height:dword
; s1p = rcx
; s2p = rdx
; dstp = r8
; stride = r9d

andMasks_SSE2 proc public frame

	width_ equ dword ptr [rbp+48]
	height equ dword ptr [rbp+56]

	push rbp
.pushreg rbp
	mov rbp,rsp
	push rbx
.pushreg rbx
	push rsi
.pushreg rsi
	push rdi
.pushreg rdi
.endprolog
	
	mov rax,rcx
	mov rbx,rdx
	mov rdx,r8
	movsxd r8,r9d
	xor rdi,rdi
	mov edi,width_
	xor rsi,rsi
	mov esi,height
	mov r10,16

yloop:
	xor rcx,rcx
xloop:
	movdqa xmm0,[rax+rcx]
	pand xmm0,[rbx+rcx]
	movdqa [rdx+rcx],xmm0
	add rcx,r10
	cmp rcx,rdi
	jl xloop

	add rax,r8
	add rbx,r8
	add rdx,r8
	dec esi
	jnz yloop

	pop rdi
	pop rsi
	pop rbx
	pop rbp

	ret
	
andMasks_SSE2 endp



;checkSceneChange_SSE2 proc s1p:dword,s2p:dword,stride:dword,width_:dword,height:dword,diffp:dword
; s1p = rcx
; s2p = rdx
; stride = r8d
; width_ = r9d

checkSceneChange_SSE2 proc public frame

	height equ dword ptr [rbp+48]
	diffp equ qword ptr [rbp+56]

	push rbp
.pushreg rbp
	mov rbp,rsp
	push rsi
.pushreg rsi
.endprolog
	
	mov rax,rcx
	mov rsi,rdx
	movsxd r8,r8d
	xor rdx,rdx
	mov edx,r9d
	xor r9,r9
	mov r9d,height
	mov r10,16
	
	pxor xmm1,xmm1
	
yloop:
	xor rcx,rcx
xloop:
	movdqa xmm0,[rax+rcx]
	psadbw xmm0,[rsi+rcx]
	paddq xmm1,xmm0

	add rcx,r10
	cmp rcx,rdx
	jl xloop

	add rax,r8
	add rsi,r8
	dec r9d
	jnz yloop


	movdqa xmm2,xmm1
	psrldq xmm1,8
	paddq xmm2,xmm1

  mov rax,diffp
	movd QWORD PTR [rax],xmm2

	pop rsi
	pop rbp

	ret
	
checkSceneChange_SSE2 endp



;VerticalBlur3_SSE2 proc srcp:dword,dstp:dword,stride:dword,width_:dword,height:dword
; srcp = rcx
; dstp = rdx
; stride = r8d
; width_ = r9d

VerticalBlur3_SSE2 proc public frame
	
	height equ dword ptr [rbp+48]

	push rbp
.pushreg rbp
	mov rbp,rsp
	push rbx
.pushreg rbx
	push rsi
.pushreg rsi
	push rdi
.pushreg rdi

	sub rsp,32
.allocstack 32
	movdqu oword ptr[rsp],xmm6
.savexmm128 xmm6,0
	movdqu oword ptr[rsp+16],xmm7
.savexmm128 xmm7,16
.endprolog

	mov rax,rcx
	mov rbx,rdx
	movsxd r8,r8d
	mov rsi,rax
	mov rdi,rax
	sub rsi,r8
	add rdi,r8
	xor rdx,rdx
	mov edx,r9d
	xor r9,r9
	mov r9d,height
	mov r10,2
	mov r11,16

	; 0x0002,for rounding
	pcmpeqb xmm6,xmm6
	psrlw xmm6,15
	psllw xmm6,1

	pxor xmm7,xmm7

	xor rcx,rcx

toploop:
	movdqa xmm0,[rax+rcx]
	pavgb xmm0,[rdi+rcx]
	movdqa [rbx+rcx],xmm0
	add rcx,r11
	cmp rcx,rdx
	jl toploop

	add rsi,r8
	add rax,r8
	add rdi,r8
	add rbx,r8

	sub r9d,r10d ; the main loop processes 2 lines fewer than the height

yloop:
	xor rcx,rcx
xloop:
	movdqa xmm0,[rsi+rcx]
	movdqa xmm1,[rax+rcx]
	movdqa xmm2,[rdi+rcx]
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
	movdqa [rbx+rcx],xmm0

	add rcx,r11
	cmp rcx,rdx
	jl xloop

	add rsi,r8
	add rax,r8
	add rdi,r8
	add rbx,r8
	dec r9d
	jnz yloop

	xor rcx,rcx

bottomloop:
	movdqa xmm0,[rsi+rcx]
	pavgb xmm0,[rax+rcx]
	movdqa [rbx+rcx],xmm0
	add rcx,r11
	cmp rcx,rdx
	jl bottomloop

	movdqu xmm7,oword ptr[rsp+16]
	movdqu xmm6,oword ptr[rsp]
	add rsp,32

	pop rdi
	pop rsi
	pop rbx
	pop rbp

	ret
	
VerticalBlur3_SSE2 endp



;HorizontalBlur3_SSE2 proc srcp:dword,dstp:dword,stride:dword,width_:dword,height:dword
; srcp = rcx
; dstp = rdx
; stride = r8d
; width_ = r9d

HorizontalBlur3_SSE2 proc public frame

	height equ dword ptr [rbp+48]

	push rbp
.pushreg rbp
	mov rbp,rsp
	sub rsp,32
.allocstack 32
	movdqu oword ptr[rsp],xmm6
.savexmm128 xmm6,0
	movdqu oword ptr[rsp+16],xmm7
.savexmm128 xmm7,16
.endprolog
	
	mov rax,rcx
	movsxd r8,r8d
	mov r10d,height
	mov r11,16

	pxor xmm7,xmm7

	pcmpeqb xmm6,xmm6
	psrlw xmm6,15
	psllw xmm6,1

yloop:
	xor rcx,rcx
xloop:
	movdqu xmm0,[rax+rcx-1]
	movdqa xmm1,[rax+rcx]
	movdqu xmm2,[rax+rcx+1]
	movdqa xmm3,xmm0
	movdqa xmm4,xmm1
	movdqa xmm5,xmm2

	punpcklbw xmm0,xmm7
	punpcklbw xmm1,xmm7
	punpcklbw xmm2,xmm7
	punpckhbw xmm3,xmm7
	punpckhbw xmm4,xmm7
	punpckhbw xmm5,xmm7

	psllw xmm1,1
	psllw xmm4,1
	paddw xmm1,xmm0
	paddw xmm4,xmm3
	paddw xmm1,xmm2
	paddw xmm4,xmm5
	psrlw xmm1,2
	psrlw xmm4,2
	packuswb xmm1,xmm4
	movdqa [rdx+rcx],xmm1

	add rcx,r11
	cmp rcx,r9
	jl xloop

	add rax,r8
	add rdx,r8
	dec r10d
	jnz yloop

	movdqu xmm7,oword ptr[rsp+16]
	movdqu xmm6,oword ptr[rsp]
	add rsp,32
	pop rbp

	ret
	
HorizontalBlur3_SSE2 endp



;HorizontalBlur6_SSE2 proc srcp:dword,dstp:dword,stride:dword,width_:dword,height:dword
; srcp = rcx
; dstp = rdx
; stride = r8d
; width_ = r9d

HorizontalBlur6_SSE2 proc public frame
	
	height equ dword ptr [rbp+48]

	push rbp
.pushreg rbp
	mov rbp,rsp
	sub rsp,112
.allocstack 112
	movdqu oword ptr[rsp],xmm6
.savexmm128 xmm6,0
	movdqu oword ptr[rsp+16],xmm7
.savexmm128 xmm7,16
	movdqu oword ptr[rsp+32],xmm8
.savexmm128 xmm8,32
	movdqu oword ptr[rsp+48],xmm9
.savexmm128 xmm9,48
	movdqu oword ptr[rsp+64],xmm10
.savexmm128 xmm10,64
	movdqu oword ptr[rsp+80],xmm11
.savexmm128 xmm11,80
	movdqu oword ptr[rsp+96],xmm12
.savexmm128 xmm12,96
.endprolog

	mov rax,rcx
	movsxd r8,r8d
	mov r10d,height
	mov r11,16

	pxor xmm12,xmm12

    ; 0x0006
    pcmpeqb xmm11,xmm11
    psrlw xmm11,14
    psllw xmm11,1

    ; 0x0008
    pcmpeqb xmm10,xmm10
    psrlw xmm10,15
    psllw xmm10,3

yloop:
	xor rcx,rcx
xloop:
	movdqu xmm0,[rax+rcx-2]
	movdqu xmm1,[rax+rcx-1]
	movdqa xmm2,[rax+rcx]
	movdqu xmm3,[rax+rcx+1]
	movdqu xmm4,[rax+rcx+2]
	movdqa xmm5,xmm0
	movdqa xmm6,xmm1
	movdqa xmm7,xmm2
	movdqa xmm8,xmm3
	movdqa xmm9,xmm4
	punpcklbw xmm0,xmm12
	punpcklbw xmm1,xmm12
	punpcklbw xmm2,xmm12
	punpcklbw xmm3,xmm12
	punpcklbw xmm4,xmm12
	punpckhbw xmm5,xmm12
	punpckhbw xmm6,xmm12
	punpckhbw xmm7,xmm12
	punpckhbw xmm8,xmm12
	punpckhbw xmm9,xmm12

	; srcp[x-2] + srcp[x+2]
	paddw xmm0,xmm4
	paddw xmm5,xmm9

	; srcp[x-1] + srcp[x+1]
	paddw xmm1,xmm3
	paddw xmm6,xmm8

	; (srcp[x-1 + srcp[x+1])*4
	psllw xmm1,2
	psllw xmm6,2

	; (srcp[x-1] + srcp[x+1])*4 + srcp[x-2] + srcp[x+2]
	paddw xmm0,xmm1
	paddw xmm5,xmm6

	; srcp[x] * 6
	pmullw xmm2,xmm11
	pmullw xmm7,xmm11
	paddw xmm0,xmm2
	paddw xmm5,xmm7

	; add 8
	paddw xmm0,xmm10
	paddw xmm5,xmm10

	; divide by 16
	psrlw xmm0,4
	psrlw xmm5,4
	packuswb xmm0,xmm5
	movdqa [rdx+rcx],xmm0

	add rcx,r11
	cmp rcx,r9
	jl xloop

	add rax,r8
	add rdx,r8
	dec r10d
	jnz yloop

	movdqu xmm12,oword ptr[rsp+96]
	movdqu xmm11,oword ptr[rsp+80]
	movdqu xmm10,oword ptr[rsp+64]
	movdqu xmm9,oword ptr[rsp+48]
	movdqu xmm8,oword ptr[rsp+32]
	movdqu xmm7,oword ptr[rsp+16]
	movdqu xmm6,oword ptr[rsp]
	add rsp,112
	pop rbp

	ret
	
HorizontalBlur6_SSE2 endp



end