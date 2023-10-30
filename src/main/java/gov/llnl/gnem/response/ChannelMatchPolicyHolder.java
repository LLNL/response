/*-
 * #%L
 * Seismic Response Processing Module
 *  LLNL-CODE-856351
 *  This work was performed under the auspices of the U.S. Department of Energy
 *  by Lawrence Livermore National Laboratory under Contract DE-AC52-07NA27344.
 * %%
 * Copyright (C) 2023 Lawrence Livermore National Laboratory
 * %%
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 * 
 *      http://www.apache.org/licenses/LICENSE-2.0
 * 
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 * #L%
 */

package gov.llnl.gnem.response;

import com.isti.jevalresp.ChannelMatchPolicy;
import com.isti.jevalresp.ChannelMatchPolicy.Policy;

/**
 *
 * @author dodge1
 */
public class ChannelMatchPolicyHolder {
    ChannelMatchPolicy policy;
    private ChannelMatchPolicyHolder() {
        policy = new ChannelMatchPolicy(Policy.FULL_MATCH);
    }
    
    public static ChannelMatchPolicyHolder getInstance() {
        return ChannelMatchPolicyHolderHolder.INSTANCE;
    }
    
    private static class ChannelMatchPolicyHolderHolder {

        private static final ChannelMatchPolicyHolder INSTANCE = new ChannelMatchPolicyHolder();
    }

    public ChannelMatchPolicy getPolicy() {
        return policy;
    }
    
    public void setPolicy(Policy policy){
        this.policy = new ChannelMatchPolicy(policy);
    }
}
