
所有点的度均为偶数的前提下，寻找欧拉回路，并把简单环拆分出来。

[ 武汉邀请赛E ](https://codeforces.com/gym/105901/problem/E)

```
void dfs(int u){
	for(int &i = cur[u]; i < t[u].size();){
		int v = t[u][i];
		if(e[HASH(u,v)] < E[HASH(u,v)]){
			e[HASH(u, v)]++, e[HASH(v, u)]++;
			i++; d[u]--, d[v]--;
			dfs(v);
		}
		else i++;
	}
	
	if(in_stk[u] == false){
		stk.push(u), in_stk[u] = true;
	}
	else{
		int pre_pos = u;
		col ++;
		while(!stk.empty()){
			int pos = stk.top();
			stk.pop();
			
			edge[HASH(pos, pre_pos)].push(col) , edge[HASH(pre_pos, pos)].push(col);
			pre_pos = pos;
			in_stk[pos] = false;
			if(pos == u) break;
		}
		
		stk.push(u), in_stk[u] = true;
	}
}
```
